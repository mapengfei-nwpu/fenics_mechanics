Python script to convert Gmsh mesh into DOLFIN format, with
support for high order elements.

Usage
-----
```bash
gmsh2dolfin [ifile].msh [-o ofile.h5]
geo2dolfin [ifile].geo [-o ofile.h5] [-t topo_dim] [-g geo_dim]
```
The output is a .h5 file containing the mesh, the entity markers
and possibly the high-order coordinates.

Remarks
-------
* It requires DOLFIN
* In principle, there is no restriction on the maximum order
  In practice, Gmsh doesn't generate 3D tetrahedral meshes
  of order > 2
* Embedding dimension is deduced from the coordinates
* Support from quad and hexa is in place but not finished

