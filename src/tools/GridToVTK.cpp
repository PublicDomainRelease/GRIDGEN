
#include "GridToVTK.h"
#include <fstream>
#include <sstream>
using namespace std;

#define ZSCALE 1

//int debug_num[4] = {0};

//reverse interpolate the z coordinates for all boxes
void reverseInterpolate(QuadTree3D* qtree);

inline void rotateOneBox(Box* b, double cx, double cy, double rotation, double* padfX, double* padfY, double* padfZ, double* padfM, double* ctrX, double* ctrY)
{
	double ldX, ldY, luX, luY, rdX, rdY, ruX, ruY, ctrx, ctry;
	b->rotate(cx, cy, rotation, ldX, ldY, luX, luY, rdX, rdY, ruX, ruY);
	b->rotateCtr(cx, cy, rotation, ctrx, ctry);

	padfX[0] = luX; padfX[1] = ruX;
	padfX[2] = rdX; padfX[3] = ldX;
	padfY[0] = luY; padfY[1] = ruY;
	padfY[2] = rdY; padfY[3] = ldY;
	ctrX[0] = ctrx; ctrY[0] = ctry;
}

/********************************************************/
//write the vtk file
void writeVTKHeader(ofstream& ofile)
{
	ofile<<"<?xml version=\"1.0\"?>\n";
	ofile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
	ofile<<"<UnstructuredGrid GhostLevel=\"0\">\n";
}
void writeVTKTail(ofstream& ofile)
{
	ofile<<"</UnstructuredGrid>\n";
	ofile<<"</VTKFile>\n";
}

void writeVTKCell(ofstream& ofile, int boxNum)
{
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	int sid = 0, id = 0;
	for(int i = 0; i < boxNum; i++)
	{
		ofile<<(sid)<<"\t"<<(sid+6)<<"\t"<<(sid+2)<<"\t"<<(sid+4)<<"\t"
		<<(sid+1)<<"\t"<<(sid+7)<<"\t"<<(sid+3)<<"\t"<<(sid+5)<<"\n";

			//ofile<<(sid+6)<<"\t"<<sid<<"\t"<<(sid+4)<<"\t"<<(sid+2)<<"\t";
			//	ofile<<(sid+7)<<"\t"<<(sid+1)<<"\t"<<(sid+5)<<"\t"<<(sid+3)<<"\n";
				sid+=8;

				id++;
				if(id%500==0)
					ofile.flush();
	}
	ofile<<"</DataArray>\n";
	//write the offset
	sid = id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < boxNum; i++)
	{
				ofile<<sid<<"\t";
				sid+=8;

				id++;
				if(id%1000==0)
					ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < boxNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%1000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";

}

template<typename Dtype>
void writeDataArray(ofstream& ofile, string type, string name, const vector<Dtype>& data)
{
	ofile<<"<DataArray type=\""<<type<<"\" Name=\""<<name<<"\" format=\"ascii\">";
	int id = 0;
	for(vector<Dtype>::const_iterator dit = data.begin(); dit != data.end(); dit++)
	{
		ofile<<(*dit)<<"\t";
		
		if(id % 1000==0)
			ofile.flush();
		id++;
	}
	ofile<<"</DataArray>\n";
}

//export modflow grid to a paraview vtk file
void writeModflowGrid2VTK(ModflowGrid* modflow, string vtkName, bool without_inactive, bool one_based_numbering)
{
	double cx, cy, rotation;
	modflow->getRotatePara(cx, cy, rotation);
	//write the points, loop over each layer
	int nlay = modflow->nlay;
	int nrow = modflow->nrow;
	int ncol = modflow->ncol;

	double padfX[4] = {0};
	double padfY[4] = {0};
	double padfZ[4] = {0};
	double padfM[4] = {0};
	double ctrX[1] = {0};
	double ctrY[1] = {0};

	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);

	/***************************************************************************************/
	int totalCell = nlay * nrow * ncol;
	ofile<<"<Piece NumberOfPoints=\""<<totalCell*8<<"\" NumberOfCells=\""<<totalCell<<"\" >";

	vector<int> nodeIds, nodeNums, layers, rows, cols,vtkids;
	vector<double> tops, bottoms, delrs, delcs;

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int id = 0, cellNum = 0;
	for(int l=0;l<nlay;l++)
	{
		for(int r=0;r<nrow;r++)
		{
			for(int c=0;c<ncol;c++)
			{
				int bid=modflow->get_nodeid(l,r,c);
				Box * box=modflow->nodegroup[bid];

				//if(box->isLeaf && box->active)
				if(box->isLeaf)
				{
					if(without_inactive && (!box->active))
						continue;

					double top= modflow->bot[ l * nrow * ncol + r * ncol + c];
					double bot= modflow->bot[(l + 1) * nrow * ncol + r * ncol + c] ;

					rotateOneBox(box, cx, cy, rotation, padfX, padfY, padfZ, padfM, ctrX, ctrY);
					for(int i = 0; i < 4; i++)
					{
						ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<top * ZSCALE<<"\t" ;		//(box->z - box->dz / 2)<<"\t";
						ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<bot * ZSCALE<<"\t";		//(box->z + box->dz / 2)<<"\t";

						id+=2;
						if(id%500==0)
							ofile.flush();
					}
					
					//store info
					if(!one_based_numbering)
					{
						nodeIds.push_back(box->id);
						nodeNums.push_back(box->number); layers.push_back(l); rows.push_back(r); cols.push_back(c);
					}
					else
					{
						nodeIds.push_back(box->id + 1);
						nodeNums.push_back(box->number + 1); layers.push_back(l + 1); rows.push_back(r + 1); cols.push_back(c + 1);
					}

					tops.push_back(top); bottoms.push_back(bot); delrs.push_back(box->dx); delcs.push_back(box->dy);
				

					vtkids.push_back(cellNum);
					cellNum++;
				}
			}
		}
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write cells
	writeVTKCell(ofile, id/8);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	writeDataArray(ofile, "Int32", "id", nodeIds);
	writeDataArray(ofile, "Int32", "nodenumber", nodeNums);
	writeDataArray(ofile, "Int32", "layer", layers);
	writeDataArray(ofile, "Int32", "row", rows);
	writeDataArray(ofile, "Int32", "col", cols);
	writeDataArray(ofile, "Int32", "vtkid", vtkids);

	writeDataArray(ofile, "Float64", "top", tops);
	writeDataArray(ofile, "Float64", "bottom", bottoms);
	writeDataArray(ofile, "Float64", "delr", delrs);
	writeDataArray(ofile, "Float64", "delc", delcs);
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();

}

void DFSWriteQDT2VTK(QuadTree3D* qtree, ofstream& ofile,  int& pntNum, Box* b, int l,int r, int c,bool one_based_numbering,double cx, double cy, double rotation,
	vector<int>* intvals, vector<double>* dbvals)
{
	char tmp[8]={'\0'};
	//if(b->isLeaf)
	if(b->isLeaf && b->active)
	{
		double padfX[4] = {0};
		double padfY[4] = {0};
		double padfZ[4] = {0};
		double padfM[4] = {0};
		double ctrX[1] = {0};
		double ctrY[1] = {0};

		rotateOneBox(b, cx, cy, rotation, padfX, padfY, padfZ, padfM, ctrX, ctrY);

		double top=qtree->top[b->number];
		double bot=qtree->bot[b->number];

		for(int i = 0; i < 4; i++)
		{
				ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<top * ZSCALE<<"\t";
				ofile<<padfX[i]<<"\t"<<padfY[i]<<"\t"<<bot * ZSCALE<<"\t";

				pntNum+=2;
				if(pntNum%500==0)
					ofile.flush();
		}
		//store info

		if(!one_based_numbering)
		{
			intvals[0].push_back(b->id); 
			intvals[1].push_back(b->number); 
			intvals[2].push_back(l); 
			intvals[3].push_back(r); 
			intvals[4].push_back(c);
		}
		else
		{
			intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
			intvals[2].push_back(l + 1); intvals[3].push_back(r + 1); intvals[4].push_back(c + 1);
		}

		dbvals[0].push_back(top); 
		dbvals[1].push_back(bot); 
		dbvals[2].push_back(b->dx); 
		dbvals[3].push_back(b->dy);

		intvals[5].push_back(pntNum);//the vtk_id
	}
	else if(b->isLeaf==false)//b is not a leaf, call this function recursively
	{
		for (int i = 0; i < 4; i++)
		{
			//if(one_based_numbering) sprintf(tmp,"%d",i+1);
			//else sprintf(tmp,"%d",i);
			DFSWriteQDT2VTK(qtree, ofile,pntNum,b->pChildren[i],l,r,c,one_based_numbering, cx, cy, rotation, intvals, dbvals);
		}
	}//end b->isLeaf
}


//write the Quadtree grid
void writeQaudtreeGrid2VTK(QuadTree3D* qtree, string vtkName)
{
	//make sure the number of the nodes are correct
	qtree->number_nodes();


	ModflowGrid * mfgrid=qtree->getModflowGrid();

	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	int totalCell = qtree->nodes();
	ofile<<"<Piece NumberOfPoints=\""<<totalCell*8<<"\" NumberOfCells=\""<<totalCell<<"\" >";

	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<double> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int pntNum = 0;
	for(int l=0;l<mfgrid->nlay;l++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(l,r,c);
				Box * box=qtree->nodegroup[bid];

				if(box->active)
					DFSWriteQDT2VTK(qtree, ofile, pntNum,  box,  l, r,  c,qtree->get_one_based_node_numbering(), cx,  cy,  rotation, intvals, dbvals);
			}
		}
	}

	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write cells
	writeVTKCell(ofile, pntNum/8);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float64", dbNames[i], dbvals[i]);
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();
}

void writeCells(ofstream& ofile,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, BoxVtx* boxVts)
{
	//ModflowGrid* mfgrid = qtree->getModflowGrid();
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	int rid = 0;
	int cellNum = allBoxes.size();
	for(vector<Box*>::iterator bit = allBoxes.begin(); bit != allBoxes.end(); ++bit)
	{
		Box* b = *bit;
		ofile<<boxVts[b->id].vertices[0]->id<<"\t"<<boxVts[b->id].vertices[1]->id<<"\t"<<boxVts[b->id].vertices[3]->id<<"\t"<<boxVts[b->id].vertices[2]->id<<"\t"
			<<boxVts[b->id].vertices[4]->id<<"\t"<<boxVts[b->id].vertices[5]->id<<"\t"<<boxVts[b->id].vertices[7]->id<<"\t"<<boxVts[b->id].vertices[6]->id<<"\n";

		rid++;
		if(rid % 5000 == 0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";

	//write the offset
	int sid, id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < cellNum; i++)
	{
		ofile<<sid<<"\t";
		sid+=8;

		id++;
		if(id%5000==0) ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < cellNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%5000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";
}
/******************************************************************************************************************
* write VTK file with sharing vertices for modflow grid
*******************************************************************************************************************/
void buildSharingVertices_mfgrid(ModflowGrid* mfgrid, vector<Vertex*>& vdata, BoxVtx* boxVts)
{
	double val;//, resx, resy;
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	int mfgrid_csize=(mfgrid->nrow) * (mfgrid->ncol);
	
	//get all the Xe and Ye
	double* Xe = mfgrid->get_local_Xe_array();
	double* Ye = mfgrid->get_local_Ye_array();

	//use the edges to locate the box corners
	int vid = 0;
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		double* zgrid = &(mfgrid->bot[mfgrid_csize * k]);
		mf_ascii_grid ag(mfgrid, zgrid, mfgrid->ncol, mfgrid->nrow);

		for(int r=0; r<=mfgrid->nrow;r++)
		{
			for(int c=0; c<=mfgrid->ncol;c++)
			{
				Point2d pos(Xe[c], Ye[r]);
				val = QuadTree3D::bilinear_interpolation(pos, &ag);

				//rotateOnePnt(pos[0], pos[1], cx, cy, rotation, resx, resy);

				Vertex* nv = new Vertex(pos[0], pos[1], val, vid++);
				vdata.push_back(nv);
			}
		}
	}

	int mfgrid_esize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
	int perrow_esize = mfgrid->ncol + 1;
	//new assign the vertex pointer to boxes
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r <	mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * box=mfgrid->nodegroup[bid];
				boxVts[box->id].vertices[0] = vdata[ k * mfgrid_esize + r * perrow_esize + c];
				boxVts[box->id].vertices[1] = vdata[ k * mfgrid_esize + (r + 1) * perrow_esize + c];
				boxVts[box->id].vertices[2] = vdata[ k * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
				boxVts[box->id].vertices[3] = vdata[ k * mfgrid_esize + r * perrow_esize + c + 1];

				boxVts[box->id].vertices[4] = vdata[ (k + 1) * mfgrid_esize + r * perrow_esize + c];
				boxVts[box->id].vertices[5] = vdata[ (k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c];
				boxVts[box->id].vertices[6] = vdata[ (k + 1) * mfgrid_esize + (r + 1) * perrow_esize + c + 1];
				boxVts[box->id].vertices[7] = vdata[ (k + 1) * mfgrid_esize + r * perrow_esize + c + 1];
			}
		}
	}

}
int collectCell_share(ModflowGrid* mfgrid,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, bool without_inactive, bool one_based_numbering)
{
	int cellNum = 0;
	for(int k=0;k<mfgrid->nlay;k++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * b=mfgrid->nodegroup[bid];
				if(b->isLeaf)
				{
					if(without_inactive && (!b->active))
						continue;

					allBoxes.push_back(b);

					//store info
					if(one_based_numbering)
					{
						intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
						intvals[2].push_back(k + 1); intvals[3].push_back(r + 1); intvals[4].push_back(c + 1);
					}
					else
					{
						intvals[0].push_back(b->id); intvals[1].push_back(b->number); 
						intvals[2].push_back(k); intvals[3].push_back(r); intvals[4].push_back(c);
					}

					double top=mfgrid->bot[k * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];
					double bot=mfgrid->bot[ (k+1) * mfgrid->nrow * mfgrid->ncol + r * mfgrid->ncol + c];

					dbvals[0].push_back((float)top); 
					dbvals[1].push_back((float)bot); 
					dbvals[2].push_back((float)(b->dx)); 
					dbvals[3].push_back((float)(b->dy));

					intvals[5].push_back(cellNum++);//the vtk_id
				}
			}
		}
	}
	return cellNum;
}

void writeModflowGrid2VTK_shareV(ModflowGrid* mfgrid, string vtkName, bool without_inactive, bool one_based_numbering)
{
	int boxNum = mfgrid->nodegroup.size();
	BoxVtx* boxVts = new BoxVtx[boxNum];
	vector<Vertex*> vdata;
	buildSharingVertices_mfgrid(mfgrid, vdata, boxVts);
		//begin to write
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	/*****************************************************************************************/
	//collect all the boxes and necessary info
	vector<Box*> allBoxes;
	
	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<float> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");

	int cellNum  = collectCell_share(mfgrid,  allBoxes, intvals, dbvals, without_inactive, one_based_numbering);

	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	ofile<<"<Piece NumberOfPoints=\""<<vdata.size()<<"\" NumberOfCells=\""<<allBoxes.size()<<"\" >";

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int flash = 0;
	for(vector<Vertex*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		Vertex* nv = *vit;
		double resx, resy;
		rotateOnePnt(nv->pos[0], nv->pos[1], cx, cy, rotation,resx, resy);
		ofile<<(float)(resx)<<"\t"<<(float)(resy)<<"\t"<<(float)(nv->pos[2] * ZSCALE)<<"\t";

		if(flash % 5000 == 0)
		{
			ofile.flush();
		}
		flash++;
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write the boxes
	writeCells(ofile,  allBoxes, intvals, dbvals, boxVts);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float32", dbNames[i], (dbvals[i]));
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();

	//clear memory
	for(vector<Vertex*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		delete (*vit);
	}
	delete[] boxVts;
}

/******************************************************************************************************************
* write VTK file with sharing vertices for quadtree
******************************************************************************************************************/
//void DFSCollectInfo(Box* b, int l, int r, int c, int& boxNum, 
void DFScollectCell_share(Box* b, int k, int r, int c, QuadTree3D* qtree, vector<Box*>& allBoxes, int& boxNum,  vector<int>* intvals, vector<float>* dbvals
	, bool one_based_numbering)
{
	if(b->isLeaf && b->active)
	{
		allBoxes.push_back(b);

		//store info
		if(one_based_numbering)
		{
			intvals[0].push_back(b->id + 1); intvals[1].push_back(b->number + 1); 
			intvals[2].push_back(k + 1); intvals[3].push_back(r + 1);intvals[4].push_back(c + 1);
		}
		else
		{
			intvals[0].push_back(b->id); intvals[1].push_back(b->number); 
			intvals[2].push_back(k); intvals[3].push_back(r);intvals[4].push_back(c);
		}

		double top=qtree->top[b->number];
		double bot=qtree->bot[b->number];

		dbvals[0].push_back((float)top); 
		dbvals[1].push_back((float)bot); 
		dbvals[2].push_back((float)(b->dx)); 
		dbvals[3].push_back((float)(b->dy));

		intvals[5].push_back(boxNum++);//the vtk_id
	}
	else if(b->isLeaf==false)
	{
		for(int i = 0; i < 4; i++)
		{
			DFScollectCell_share(b->pChildren[i], k, r, c, qtree, allBoxes, boxNum, intvals, dbvals, one_based_numbering);
		}
	}
}
int collectCell_share(QuadTree3D* qtree,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, bool one_based_numbering)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	int cellNum = 0;
	for(int k=0;k<mfgrid->nlay;k++)
	{
		for(int r=0;r<mfgrid->nrow;r++)
		{
			for(int c=0;c<mfgrid->ncol;c++)
			{
				int bid=mfgrid->get_nodeid(k,r,c);
				Box * box=qtree->nodegroup[bid];
				if(box->active)
					DFScollectCell_share(box, k, r, c, qtree, allBoxes, cellNum, intvals, dbvals, one_based_numbering);
			}
		}
	}
	return cellNum;
}

Vertex* getMid(Vertex* v1, Vertex* v2, vector<Vertex*>& vdata, int& vid)
{
	Vertex* nv = new Vertex( (v1->pos[0] + v2->pos[0])/2, (v1->pos[1] +v2->pos[1])/2, (v1->pos[2] + v2->pos[2])/2, vid++);
	vdata.push_back(nv);

	return nv;
}

void rotateOnePnt(double x, double y, double cx, double cy, double rotation,double& resx, double& resy)
{
	resx= cx +  (x - cx) * cos(rotation) - (y - cy)* sin(rotation);
	resy = cy + (y - cy)*cos(rotation) + (x - cx)*sin(rotation);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//fit a plane using 3 points
double fitZ(double** vv, double x, double y)
{
	double a = (vv[1][1] - vv[0][1]) * (vv[2][2] - vv[0][2]) - (vv[2][1] - vv[0][1]) * (vv[1][2] - vv[0][2]);
	double b = (vv[1][2] - vv[0][2]) * (vv[2][0] - vv[0][0]) - (vv[2][2] - vv[0][2]) * (vv[1][0] - vv[0][0]);
	double c = (vv[1][0] - vv[0][0]) * (vv[2][1] - vv[0][1]) - (vv[2][0] - vv[0][0]) * (vv[1][1] - vv[0][1]);
	double d = -(a * vv[0][0] + b * vv[0][1] + c * vv[0][2]);
	
	assert(abs(c)>=1e-10);
	double z = - (a * x + b * y + d) / c;
	
	return z;
}


///////////////////////////////////////////////////////////////////////////////////////////
struct Share_E;
struct Share_Box;

/*
	edge order
		3
		|
		|
 0	---------	2
		|
		|
		1
*/
struct Share_V
{
	Share_V()
	{
		for(int i = 0; i < 4; i++)
		{
			e[i] = NULL; ctr_e[i] = NULL;
		}
		for(int i = 0; i < 8; i++)
			box[i] = NULL;
		layer = num = -1;
		is_midPnt = false;
		parent_e = NULL;
	}
	Share_V(double x, double y)
	{
		for(int i = 0; i < 4; i++)
		{
			e[i] = NULL;
		}
		for(int i = 0; i < 8; i++)
			box[i] = NULL;
		layer = num = -1;
		pos[0] = x; pos[1] = y;
		is_midPnt = false;
		parent_e = NULL;
	}

	Point3d pos;
	Share_E* e[4];
	Box* box[8];
	int layer, num;
	Share_E* ctr_e[4];//only valid for splitted box
	bool is_midPnt;
	Share_E* parent_e;//only valid when is_midPnt == true
};
/*
 vertex order: alway from left to right, top to bottom
 0 --------- 1

 0
 |
 |
 |
 1

*/
struct Share_E
{
	Share_E(Share_V* v1, Share_V* v2)
	{
		parent = NULL;
		child[0] = child[1] = NULL;
		v[0] = v1; v[1] = v2;
		box[0] = box[1] = box[2] = box[3] = NULL;
		depth = 0;
	}
	Share_E* parent;
	Share_E* child[2];
	Share_V* v[2]; //the two attached vertices
	Box* box[4];//the two attached boxes, 0 for left or up, 1 for right or down(by indxe from 0 to nrow)
	int depth;
	//the box center enclosed by these two edges
	map<Share_E*, Share_V*> ee2v;
};

/*
*	edges
		3
	____________
	|			|
0	|			|	2
	|			|
	------------
		1
*/
struct Share_Box
{
	Share_Box()
	{
		for(int i = 0; i < 8; i++)
		{
			e[i] = NULL; //box[i] = NULL; 
		}
	}
	Share_E* e[8];
	//Box* box[8]; //assuming that the depth difference between two nearby boxes is at most 1
};
void assignChildBoxEdge(Box* box, Share_E* es[8], Share_Box* boxdata)
{
	for(int i = 0; i < 2; i++)
	{
		Share_Box& s_box = boxdata[box->pChildren[0]->id];

		s_box.e[i * 4] = boxdata[box->id].e[i * 4]->child[0];
		s_box.e[i * 4 + 1] = es[i * 4];
		s_box.e[i * 4 + 2] = es[i * 4 + 3];
		s_box.e[i * 4 + 3] = boxdata[box->id].e[i * 4 + 3]->child[0];
	}
	for(int i = 0; i < 2; i++)
	{
		Share_Box& s_box = boxdata[box->pChildren[3]->id];

		s_box.e[i * 4] = boxdata[box->id].e[i * 4]->child[1];
		s_box.e[i * 4 + 1] = boxdata[box->id].e[i * 4 + 1]->child[0];
		s_box.e[i * 4 + 2] = es[i * 4 + 1];
		s_box.e[i * 4 + 3] = es[i * 4];
	}
	for(int i = 0; i < 2; i++)
	{
		Share_Box& s_box = boxdata[box->pChildren[2]->id];

		s_box.e[i * 4] = es[i * 4 + 1];
		s_box.e[i * 4 + 1] = boxdata[box->id].e[i * 4 + 1]->child[1];
		s_box.e[i * 4 + 2] = boxdata[box->id].e[i * 4 + 2]->child[1];
		s_box.e[i * 4 + 3] = es[i * 4 + 2];
	}
	for(int i = 0; i < 2; i++)
	{
		Share_Box& s_box = boxdata[box->pChildren[1]->id];

		s_box.e[i * 4] = es[i * 4 + 3];
		s_box.e[i * 4 + 1] = es[i * 4 + 2];
		s_box.e[i * 4 + 2] = boxdata[box->id].e[i * 4 + 2]->child[0];
		s_box.e[i * 4 + 3] = boxdata[box->id].e[i * 4 + 3]->child[1];
	}

}
//use DFS to build vertices
void DFSBuildShareV(Box* box, vector<Share_V*>& vdata, vector<Share_E*>& edata, Share_Box* boxdata)
{
	if(box->isLeaf) return;
	//break the edge into two pieces
	for(int k = 0; k < 2; k++)
	{
		for(int t = 0; t < 4; t++)
		{
			Share_E* te = boxdata[box->id].e[k * 4 + t];
			if(te->child[0] == NULL)
			{
				//break this edge
				Share_V* nv = new Share_V((te->v[0]->pos[0] + te->v[1]->pos[0])/2, (te->v[0]->pos[1] + te->v[1]->pos[1])/2);
				nv->layer = te->v[0]->layer;
				vdata.push_back(nv);

				Share_E* ne1 = new Share_E(te->v[0], nv);
				Share_E* ne2 = new Share_E(nv, te->v[1]);
				edata.push_back(ne1); edata.push_back(ne2);
				te->child[0] = ne1; te->child[1] = ne2;
				ne1->parent = te; ne2->parent = te;

				//update the edges for v
				if(t == 0 || t == 2) 
				{
					te->v[0]->e[1] = ne1; te->v[1]->e[3] = ne2;
					nv->e[3] = ne1; nv->e[1] = ne2;
				}
				else
				{
					te->v[0]->e[2] = ne1; te->v[1]->e[0] = ne2;
					nv->e[0] = ne1; nv->e[2] = ne2;
				}
			}
		}
	}
	//construct the center
	Share_V* vs[2] = {NULL};
	Share_E* es[8] = {NULL};

	for(int k = 0; k < 2; k++)
	{
		if(boxdata[box->id].e[k * 4]->ee2v.find(boxdata[box->id].e[k * 4 + 2]) != boxdata[box->id].e[k * 4]->ee2v.end())
		{
			vs[k] = boxdata[box->id].e[k * 4]->ee2v[boxdata[box->id].e[k * 4 + 2]];
			for(int i = 0; i < 4; i++) es[k * 4 + i] = vs[k * 4 + i]->ctr_e[i];
		}
		else
		{
			vs[k] = new Share_V(box->x, box->y);
			vdata.push_back(vs[k]);
			vs[k]->layer = boxdata[box->id].e[k * 4]->v[0]->layer;
			for(int i = 0; i < 4; i++)
			{
				if(i == 0 || i == 3)
					es[k * 4 + i] = new Share_E(boxdata[box->id].e[k * 4 + i]->child[0]->v[1], vs[k]);
				else
					es[k * 4 + i] = new Share_E(vs[k], boxdata[box->id].e[k * 4 + i]->child[0]->v[1]);
				vs[k]->ctr_e[i] = es[k * 4 + i];
			}		
			es[k * 4]->ee2v.insert(make_pair(es[k * 4 + 2], vs[k])); es[k * 4 + 2]->ee2v.insert(make_pair(es[k * 4], vs[k]));
			es[k * 4 + 1]->ee2v.insert(make_pair(es[k * 4 + 3], vs[k])); es[k * 4 + 3]->ee2v.insert(make_pair(es[k * 4 + 1], vs[k]));

			//update
			for(int i = 0; i < 4; i++)
			{
				vs[k]->e[i] = es[k * 4 + i];
			}
			es[k * 4]->v[0]->e[2] = es[k * 4];
			es[k * 4 + 1]->v[1]->e[3] = es[k * 4 + 1];
			es[k * 4 + 2]->v[1]->e[0] = es[k * 4 + 2];
			es[k * 4 + 3]->v[0]->e[1] = es[k * 4 + 3];
		}
	}

	//assign the edges to child boxes
	assignChildBoxEdge(box, es, boxdata);

	for(int i = 0; i < 4; i++)
	{
		DFSBuildShareV(box->pChildren[i], vdata, edata, boxdata);
	}
}
//determine: edges for v, boxes for e, nearby boxes for box
void determinVEB(QuadTree3D* qtree, vector<Share_V*>& vdata, vector<Share_E*>& edata, Share_Box* boxdata)
{
	int bsize = qtree->nodegroup.size();
	for(int i = 0; i < bsize; i++)
	{
		Box* b = qtree->nodegroup[i];
		if(!(b->isLeaf && b->active))
			continue;

		Share_Box& s_b = boxdata[b->id];
		//determine boxes for e
		for(int k = 0; k < 2; k++)
		{
			for(int t = 0; t < 4; t++)
			{
				if(t == 0 || t == 3)
					s_b.e[k * 4 + t]->box[k == 0 ? 3 : 1] = b;
				else
				{
					s_b.e[k * 4 + t]->box[k == 0 ? 2 : 0] = b;
				}
			}
		}
		 
		//determine boxes for v
		for(int k = 0; k < 2; k++)
		{
			s_b.e[k * 4 + 0]->v[0]->box[(1-k) * 4 + 2] = b;
			if(s_b.e[k * 4]->child[0] != NULL)
				s_b.e[k * 4]->child[0]->v[1]->box[(1-k) * 4 + 2] = b;

			s_b.e[k * 4 + 1]->v[0]->box[(1-k) * 4 + 3] = b;
			if(s_b.e[k * 4 + 1]->child[0] != NULL)
				s_b.e[k * 4 + 1]->child[0]->v[1]->box[(1-k) * 4 + 3] = b;

			s_b.e[k * 4 + 2]->v[1]->box[(1-k) * 4 + 0] = b;
			if(s_b.e[k * 4 + 2]->child[0] != NULL)
				s_b.e[k * 4 + 2]->child[0]->v[1]->box[(1-k) * 4 + 0] = b;

			s_b.e[k * 4 + 3]->v[1]->box[(1-k) * 4 + 1] = b;
			if(s_b.e[k * 4 + 3]->child[0] != NULL)
				s_b.e[k * 4 + 3]->child[0]->v[1]->box[(1-k) * 4 + 1] = b;
		}

		//determine the mid point 
		for(int i = 0; i < 8; i++)
		{
			if(s_b.e[i]->child[0] != NULL)
			{
				s_b.e[i]->child[0]->v[1]->is_midPnt = true;
				s_b.e[i]->child[0]->v[1]->parent_e = s_b.e[i];
			}
		}
	}

	//update 
	for(vector<Share_E*>::iterator eit = edata.begin(); eit != edata.end(); ++eit)
	{
		Share_E* te = *eit;
		for(int i = 0; i < 4; i++)
		{
			if(te->box[i] == NULL && te->parent != NULL && te->parent->box[i] != NULL)
				te->box[i] = te->parent->box[i];
		}
	}
}


//build the vertices data here
void buildSharingVertices(QuadTree3D* qtree, vector<Share_V*>& vdata, vector<Share_E*>& edata, Share_Box* boxdata)
{
	int bsize = qtree->nodegroup.size();

	ModflowGrid* mfgrid = qtree->getModflowGrid();
	double* Xe = mfgrid->get_local_Xe_array();
	double* Ye = mfgrid->get_local_Ye_array();

	//////////////////////////////////////////////////////////////////////////////
	//initialize sharing v
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		for(int r = 0; r <= mfgrid->nrow; r++)
		{
			for(int c = 0; c <= mfgrid->ncol; c++)
			{
				Share_V* nv = new Share_V();
				nv->pos[0] = Xe[c]; nv->pos[1] = Ye[r];
				nv->layer = k;
				vdata.push_back(nv);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	//initializing the sharing edge
	int mfgrid_vsize = (mfgrid->nrow + 1) * (mfgrid->ncol + 1);
	int perrow_vsize = mfgrid->ncol + 1;
	for(int k = 0; k <= mfgrid->nlay; k++)
	{
		//initializing the row edges
		for(int r = 0; r <= mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				Share_V* v1 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c];
				Share_V* v2 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c + 1];
				Share_E* e = new Share_E(v1, v2);
				v1->e[2] = e; v2->e[0] = e;
				edata.push_back(e);
			}
		}

		//initializing the col edges
		for(int c = 0; c <= mfgrid->ncol; c++)
		{
			for(int r = 0; r < mfgrid->nrow; r++)
			{
				Share_V* v1 = vdata[ k * mfgrid_vsize + r * perrow_vsize + c];
				Share_V* v2 = vdata[ k * mfgrid_vsize + (r + 1) * perrow_vsize + c];
				Share_E* e = new Share_E(v1, v2);
				v1->e[1] = e; v1->e[3] = e;
				edata.push_back(e);
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	//initialize the root box
	int mfgrid_esize = mfgrid->ncol * ( mfgrid->nrow + 1) + mfgrid->nrow * (mfgrid->ncol + 1);
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r < mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid = mfgrid->get_nodeid(k, r, c);
				Box* box = qtree->nodegroup[bid];
				
				for(int t = 0; t < 2; t++)
				{
					int baseid = (k + t) * mfgrid_esize;
					int id1 = baseid + (mfgrid->nrow + 1) * mfgrid->ncol + c * mfgrid->nrow + r;
					int id2 = baseid + (r + 1) * mfgrid->ncol + c;
					int id3 = baseid + (mfgrid->nrow + 1) * mfgrid->ncol + (c + 1) * mfgrid->nrow + r;
					int id4 = baseid + r * mfgrid->ncol + c;
					
					boxdata[box->id].e[t * 4] = edata[id1];
					boxdata[box->id].e[t * 4 + 1] = edata[id2];
					boxdata[box->id].e[t * 4 + 2] = edata[id3];
					boxdata[box->id].e[t * 4 + 3] = edata[id4];
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////
	//use DFS to build new vertices, edges and boxes
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int r = 0; r < mfgrid->nrow; r++)
		{
			for(int c = 0; c < mfgrid->ncol; c++)
			{
				int bid = mfgrid->get_nodeid(k, r, c);
				Box* box = qtree->nodegroup[bid];

				DFSBuildShareV(box, vdata, edata, boxdata);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//determine: edges for v, boxes for e
	determinVEB(qtree, vdata, edata, boxdata);
}

int number_Share_V(QuadTree3D* qtree, vector<Share_V*>& vdata, Share_Box* boxdata, vector<Share_V*>& node_vs)
{
	int bsize = qtree->nodegroup.size();
	////////////////////////////////////////////////////////////////////////////// 
	//get all the vertices that will be shared by leaf & active boxes
	int num = 0;
	for(int i = 0; i < bsize; i++)
	{
		Box* b = qtree->nodegroup[i];
		if(!(b->isLeaf && b->active))
			continue;
		Share_Box& s_b = boxdata[b->id];
		for(int k = 0; k < 8; k++)
		{
			for(int t = 0; t <2; t++)
			{
				if(s_b.e[k]->v[t]->num == -1)
				{
					node_vs.push_back(s_b.e[k]->v[t]);
					s_b.e[k]->v[t]->num = num++;
				}
			}
		}
	}
	return num;
}

int getDepth(int n[4])
{
	if(n[0] == n[1] && n[1] == n[2] && n[2] == n[3])
		return -1;
	for(int i = 0; i < 4; i++)
	{
		if(n[i] > n[(i+1)%4])
			return i;
	}
	assert(false);
	return -1;
}
bool no_collinear_add(set<Box*>& bset, Box* b)
{
	if(bset.size() == 2)
	{
		set<Box*>::iterator bit = bset.begin();

		Box* b1 = *(bset.begin());
		Box* b2 = *(++bit);
		if(bset.find(b) == bset.end())
		{
			double co_eff = abs((b2->x - b1->x) * (b->y - b1->y) - (b->x - b1->x) * (b2->y - b1->y));
			if( co_eff < 1e-6)
				return false;
		}
	}
	bset.insert(b);
	return true;
}
void getVSet(Share_Box& s_box, Share_V* v[4], bool isTop)
{
	int prefix = (isTop ? 0 : 4);
		v[0] = s_box.e[prefix + 0]->v[0]; v[1] = s_box.e[prefix + 1]->v[0]; v[2] = s_box.e[prefix + 2]->v[1]; v[3] = s_box.e[prefix + 3]->v[1];
}
void getBSet3(set<Box*>& bset, Share_Box* boxdata, bool isTop)
{
	while(bset.size() < 3)
	{
		for(set<Box*>::iterator bit = bset.begin(); bit != bset.end(); ++bit)
		{
			Box* tb = *bit;
			Share_Box& s_box = boxdata[tb->id];
			Share_V* v[4]={NULL};

			getVSet(s_box, v, isTop);
			for(int i = 0; i < 4; i++)
			{
				Share_V* stv = v[i];
				for(int j = 0; j < 4; j++)
				{
					Box* b = NULL;
					if(isTop)
						b = stv->box[4 + j];
					else
						b = stv->box[j];
					if(b != NULL && b->isLeaf && b->active)
					{
						bool tmp_b = no_collinear_add(bset, b);
						if(bset.size() == 3)
							return;
					}
				}
			}
		}
	}
}

//interpolate
double interpolate_zcoord(QuadTree3D* qtree, Share_V* v, Share_Box* boxdata, bool isTop)
{
	double* zcoord = (isTop ? qtree->top : qtree->bot);
	int prefix_bid = (isTop ? 2 : 0);
	int prefix_eid = 0;//(isTop ? 0 : 4);

	set<Box*> difset;
	Box* s[4] = {NULL};
	Box* attach_b[4] = {NULL};

	if(isTop)
	{
		for(int i = 0; i < 4; i++)
			s[i] = v->box[4 + i];
	}
	else
	{
		for(int i = 0; i < 4; i++)
			s[i] = v->box[i];
	}

	double val = 0.0;
	int num = 0;
	for(int i = 0; i < 4; i++)
		if(s[i] != NULL && s[i]->isLeaf && s[i]->active)
		{
			if(difset.find(s[i]) == difset.end())
			{
				difset.insert(s[i]);
				attach_b[num] = s[i];
				num++;
			}
		}

	//debug_num[num - 1]++;
	//if(num == 2||num == 3)
	//{
	//	cout<<"debugging...\n";
	//}

	if(num == 4)
	{
		int depth[4];
		for(int i = 0; i < 4; i++)
			depth[i] = attach_b[i]->depth;

		 int did = getDepth(depth);

		 //all 4 boxes are with same depth
		 if(did == -1)
		 {
			for(int i = 0; i < 4; i++)
				val += zcoord[attach_b[i]->number];
			return val/4.0;
		 }
		
		 //1 is smaller than other 3
		 /*
			|	   |
			|-------------
		0	|  did |
	----------------------
		1	|	2
			|

		 */
		 double v1 = (2 * zcoord[attach_b[did]->number] + zcoord[attach_b[(did+1)%4]->number])/3;
		 double v2 = (zcoord[attach_b[(did + 2)%4]->number] + zcoord[attach_b[(did+3)%4]->number])/2;
		 val = (2 * v1 + v2)/3;
	}
	else if(num == 3)
	{
		double** vv = new double*[3];
		for(int i = 0; i < 3; i++)
		{
			vv[i] = new double[3];
			memset(vv[i], 0, sizeof(double) * 3);
		}

		val = 0;
		int tid = 0;
		for(int i = 0;  i < 4; i++)
		{
			if(attach_b[i] != NULL && attach_b[i]->isLeaf && attach_b[i]->active)
			{
				vv[tid][0] = attach_b[i]->x, vv[tid][1] = attach_b[i]->y; vv[tid][2] = zcoord[attach_b[i]->number];
				tid++;
			}
		}

		val = fitZ(vv, v->pos[0], v->pos[1]);
		for(int i = 0; i < 3; i++)
			delete[] vv[i];
		delete[] vv;
	}
	else if(num == 2)
	{
		double len1 = (v->pos[0] - attach_b[0]->x) * ( v->pos[0] - attach_b[0]->x) + (v->pos[1] - attach_b[0]->y) * (v->pos[1] - attach_b[0]->y);
		double len2 = (v->pos[0] - attach_b[1]->x) * ( v->pos[0] - attach_b[1]->x) + (v->pos[1] - attach_b[1]->y) * (v->pos[1] - attach_b[1]->y);
		val = (len1 * zcoord[attach_b[1]->number] + len2 * zcoord[attach_b[0]->number]) / (len1 + len2); 
	}
	else
	{
		val = zcoord[attach_b[0]->number];
	}

	return val;
}

//do interpolation for vertices using attached boxes
void interpolate_shareV(QuadTree3D* qtree, Share_Box* boxdata, vector<Share_V*>& node_vs)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	for(vector<Share_V*>::iterator vit = node_vs.begin(); vit != node_vs.end(); ++vit)
	{
		Share_V* tv = *vit;
		if(tv->is_midPnt)
			continue;

		//interpolate from the top
		if(tv->layer == 0)
		{
			tv->pos[2] = interpolate_zcoord(qtree, tv, boxdata, true);
		}
		else if(tv->layer == mfgrid->nlay)
		{
			tv->pos[2] = interpolate_zcoord(qtree, tv, boxdata, false);
		}
		else
		{
			tv->pos[2] = (interpolate_zcoord(qtree, tv, boxdata, true) + interpolate_zcoord(qtree, tv, boxdata, false))/2;
		}
	}

	//now we process the midPnt
	for(vector<Share_V*>::iterator vit = node_vs.begin(); vit != node_vs.end(); ++vit)
	{
		Share_V* tv = *vit;
		if(!(tv->is_midPnt))
			continue;

		tv->pos[2] = (tv->parent_e->v[0]->pos[2] + tv->parent_e->v[1]->pos[2]) / 2;
	}

	//print debug info
	//cout<<"debug info in interpolation: \n";
	//for(int i = 0; i < 4; i++)
	//{
	//	cout<<debug_num[i]<<"\n";
	//}
	//cout<<endl;

}

//build the mfgrid ascii grid
void build_mf_ascii(QuadTree3D* qtree, vector<mf_ascii_grid*>& top_ags, vector<mf_ascii_grid*>& bot_ags)
{
	ModflowGrid* mfgrid = qtree->getModflowGrid();

	int csize = mfgrid->nrow * mfgrid->ncol;
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		double* top_zgrid = &(mfgrid->top[csize * k]);
		double* bot_zgrid = &(mfgrid->bot[csize * k]);
		mf_ascii_grid* top_ag = new mf_ascii_grid(mfgrid, top_zgrid, mfgrid->ncol, mfgrid->nrow); 
		mf_ascii_grid* bot_ag = new mf_ascii_grid(mfgrid, bot_zgrid, mfgrid->ncol, mfgrid->nrow);
		top_ags.push_back(top_ag);
		bot_ags.push_back(bot_ag);
	}
}

void writeCells(ofstream& ofile,  vector<Box*>& allBoxes, vector<int>* intvals, vector<float>* dbvals, Share_Box* boxdata)
{
	//ModflowGrid* mfgrid = qtree->getModflowGrid();
	/***************************************************************************************/
	//write the cells
	ofile<<"<Cells>\n";
	//write the connectivity
	ofile<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	int rid = 0;
	Share_V* vids[8];
	int cellNum = allBoxes.size();
	for(vector<Box*>::iterator bit = allBoxes.begin(); bit != allBoxes.end(); ++bit)
	{
		Box* b = *bit;
		Share_Box& s_box = boxdata[b->id];
		vids[0] = s_box.e[0]->v[0]; vids[1] = s_box.e[0]->v[1]; vids[2] = s_box.e[1]->v[1]; vids[3] = s_box.e[2]->v[0];
		vids[4] = s_box.e[4]->v[0]; vids[5] = s_box.e[4]->v[1]; vids[6] = s_box.e[5]->v[1]; vids[7] = s_box.e[6]->v[0];

		ofile<<vids[0]->num<<"\t"<<vids[1]->num<<"\t"<<vids[3]->num<<"\t"<<vids[2]->num<<"\t"
			<<vids[4]->num<<"\t"<<vids[5]->num<<"\t"<<vids[7]->num<<"\t"<<vids[6]->num<<"\n";

		rid++;
		if(rid % 5000 == 0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";

	//write the offset
	int sid, id = 0;
	sid = 8;
	ofile<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(int i =0; i < cellNum; i++)
	{
		ofile<<sid<<"\t";
		sid+=8;

		id++;
		if(id%5000==0) ofile.flush();
	}
	ofile<<"\n";
	ofile<<"</DataArray>\n";
	//write the types
	id = 0;
	ofile<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for(int i = 0; i < cellNum; i++)
	{
		ofile<<"11\t";
		id++;
		if(id%5000==0)
			ofile.flush();
	}
	ofile<<"</DataArray>\n";
	ofile<<"</Cells>\n";
}

//write the Quadtree grid with shared vertices among neighboring cells
void writeQaudtreeGrid2VTK_shareV(QuadTree3D* qtree, string vtkName)
{
	//make sure the number of the nodes are correct
	qtree->number_nodes();

	int boxNum = qtree->nodegroup.size();

	vector<Share_V*> vdata;
	vector<Share_E*> edata;
	Share_Box* boxdata = new Share_Box[boxNum];

	//build the vertices data here
	buildSharingVertices(qtree, vdata, edata, boxdata);

	//numbering vertices
	vector<Share_V*> node_vertices;
	int share_vnum = number_Share_V(qtree, vdata, boxdata, node_vertices);

	//first reversely interpolate the modflow grid top and bot
	//reverseInterpolate(qtree);
	
	//interpret the z coordinates
	interpolate_shareV(qtree,  boxdata, node_vertices); //node_vertices);

	
	ModflowGrid* mfgrid = qtree->getModflowGrid();
	//begin to write
	double cx, cy, rotation;
	mfgrid->getRotatePara(cx, cy, rotation);

	/*****************************************************************************************/
	//collect all the boxes and necessary info
	vector<Box*> allBoxes;
	
	vector<string> intNames, dbNames;
	vector<int> intvals[6];// nodeid, nodeNums, layers, rows, cols, vtkid;
	intNames.push_back("nodeid");
	intNames.push_back("nodenumber"); 
	intNames.push_back("layer"); 
	intNames.push_back("row"); 
	intNames.push_back("col");
	intNames.push_back("vtkid");

	vector<float> dbvals[4];//tops, bottoms, delrs, delcs;
	dbNames.push_back("top"); 
	dbNames.push_back("bottom"); 
	dbNames.push_back("delr"); 
	dbNames.push_back("delc");

	int cellNum  = collectCell_share(qtree,  allBoxes, intvals, dbvals, qtree->get_one_based_node_numbering());

	//write the points, loop over each layer
	stringstream ss0;
	ss0<<vtkName.c_str();
	ss0<<".vtu";
	ofstream ofile(ss0.str().c_str());
	writeVTKHeader(ofile);
	
	/***************************************************************************************/
	ofile<<"<Piece NumberOfPoints=\""<<node_vertices.size()<<"\" NumberOfCells=\""<<allBoxes.size()<<"\" >";

	//write the points
	ofile<<"<Points>\n";
	ofile<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	int flash = 0;
	for(vector<Share_V*>::iterator vit = node_vertices.begin(); vit != node_vertices.end(); ++vit)
	{
		Share_V* sv = *vit;
		double resx, resy;
		rotateOnePnt(sv->pos[0], sv->pos[1], cx, cy, rotation,resx, resy);

		ofile<<(float)(resx)<<"\t"<<(float)(resy)<<"\t"<<(float)(sv->pos[2] * ZSCALE)<<"\t";

		if(flash % 5000 == 0)
		{
			ofile.flush();
		}
		flash++;
	}

	ofile<<"\n";
	ofile<<"</DataArray>\n";
	ofile<<"</Points>\n";

	//write the boxes
	writeCells(ofile,  allBoxes, intvals, dbvals, boxdata);

	//write the attributes
	ofile<<"<CellData Scalars=\"scalars\">\n";
	int intSize = intNames.size();
	for(int i = 0; i < intSize; i++)
	{
		writeDataArray(ofile, "Int32", intNames[i], intvals[i]);
	}

	int dbSize = dbNames.size();
	for(int i = 0; i < 4; i++)
	{
		writeDataArray(ofile, "Float32", dbNames[i], (dbvals[i]));
	}
	ofile<<"</CellData>\n";

	ofile<<"</Piece>\n";
	writeVTKTail(ofile);

	ofile.close();


	//clear memory
	for(vector<Share_V*>::iterator vit = vdata.begin(); vit != vdata.end(); ++vit)
	{
		delete (*vit);
	}
	for(vector<Share_E*>::iterator eit = edata.begin(); eit != edata.end(); ++eit)
	{
		delete (*eit);
	}
	delete[] boxdata; 
}
 

//DFS add
void DFSAdd(Box* b, list<Box*>& blist)
{
    blist.push_back(b);

	if(b->isLeaf)
		return;

	for(int i = 0; i < 4; i++)
	{
		DFSAdd(b->pChildren[i], blist);
	}
}

//reverse interpolate the z coordinates for all boxes
void reverseInterpolate(QuadTree3D* qtree)
{
	int bid, bsize = qtree->nodegroup.size();
	list<Box*> blist;
	Box* b;

	double* top = new double[bsize];
	double* bot = new double[bsize];
	memset(top, 0, sizeof(double) * bsize);
	memset(bot, 0, sizeof(double) * bsize);

	ModflowGrid* mfgrid = qtree->getModflowGrid();
	//create a stack

	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int j = 0; j < mfgrid->nrow; j++)
		{
			for(int i = 0; i < mfgrid->ncol; i++)
			{
				bid = mfgrid->get_nodeid(k, j, i); 

				b = qtree->nodegroup[bid];
				DFSAdd(b, blist);				
			}
		}
	}

	//recursively use children to interpolate parent
	for(list<Box*>::reverse_iterator bit = blist.rbegin(); bit != blist.rend(); ++bit)
	{
		b = *bit;
		if(b->isLeaf)
		{
			top[b->id] = qtree->top[b->number];
			bot[b->id] = qtree->bot[b->number];
			continue;
		}
		//interpolate
		top[b->id] = (top[b->pChildren[0]->id] + top[b->pChildren[1]->id] + top[b->pChildren[2]->id] + top[b->pChildren[3]->id]) / 4.0;
		bot[b->id] = (bot[b->pChildren[0]->id] + bot[b->pChildren[1]->id] + bot[b->pChildren[2]->id] + bot[b->pChildren[3]->id]) / 4.0;
	}

	//copy back to modflow grid top and bottom
	for(int k = 0; k < mfgrid->nlay; k++)
	{
		for(int j = 0; j < mfgrid->nrow; j++)
		{
			for(int i = 0; i < mfgrid->ncol; i++)
			{
				bid = mfgrid->get_nodeid(k, j, i); 
				b = qtree->nodegroup[bid];
				mfgrid->top[bid] = top[b->id];
				mfgrid->bot[bid] = bot[b->id];
			}
		}
	}

	//clear memory
	 delete[] top;
	 delete[] bot;

}