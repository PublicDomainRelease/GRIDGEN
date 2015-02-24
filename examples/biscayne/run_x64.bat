rem build the quadtree grid grid02qtg
..\..\bin\gridgen_x64.exe buildqtg action01_buildqtg.dfn

rem write modflow-usg data for both grids
..\..\bin\gridgen_x64.exe grid02qtg-to-usgdata action02_writeusgdata.dfn

rem create shapefiles for both grids
..\..\bin\gridgen_x64.exe grid01mfg-to-polyshapefile action03_shapefile.dfn
..\..\bin\gridgen_x64.exe grid02qtg-to-polyshapefile action03_shapefile.dfn
..\..\bin\gridgen_x64.exe grid01mfg-to-pointshapefile action03_shapefile.dfn
..\..\bin\gridgen_x64.exe grid02qtg-to-pointshapefile action03_shapefile.dfn

rem intersect features with quadtree grid
..\..\bin\gridgen_x64.exe canal_grid02qtg_lay1_intersect action04_intersect.dfn
..\..\bin\gridgen_x64.exe poly_grid02qtg_lay1_intersect action04_intersect.dfn
..\..\bin\gridgen_x64.exe chd_grid02qtg_lay1_intersect action04_intersect.dfn
..\..\bin\gridgen_x64.exe lu2008_grid02qtg_lay1_intersect action04_intersect.dfn

rem create a vtk file of the grid
..\..\bin\gridgen_x64.exe grid01mfg-to-vtkfile action05_vtkfile.dfn
..\..\bin\gridgen_x64.exe grid02qtg-to-vtkfile action05_vtkfile.dfn


pause
