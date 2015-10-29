rem build the quadtree grid grid02qtg
..\..\bin\gridgen.exe buildqtg action01_buildqtg.dfn

rem write modflow-usg data
..\..\bin\gridgen.exe grid02qtg-to-usgdata action02_writeusgdata.dfn

rem create shapefiles for both grids
..\..\bin\gridgen.exe grid01mfg-to-polyshapefile action03_shapefile.dfn
..\..\bin\gridgen.exe grid02qtg-to-polyshapefile action03_shapefile.dfn
..\..\bin\gridgen.exe grid01mfg-to-pointshapefile action03_shapefile.dfn
..\..\bin\gridgen.exe grid02qtg-to-pointshapefile action03_shapefile.dfn

rem intersect features with quadtree grid
..\..\bin\gridgen.exe canal_grid02qtg_lay1_intersect action04_intersect.dfn
..\..\bin\gridgen.exe chd_grid02qtg_lay1_intersect action04_intersect.dfn

rem create a vtk file of the grid
..\..\bin\gridgen.exe grid01mfg-to-vtkfile action05_vtkfile.dfn
..\..\bin\gridgen.exe grid02qtg-to-vtkfile action05_vtkfile.dfn
..\..\bin\gridgen.exe grid02qtg-to-vtkfilesv action05_vtkfile.dfn

pause
