// =====
// specs
// -----
dfocal  = 3.7; // [cm]
l_endo  = 0.4;
l_epi   = 0.7;
mu_base = 120.0;

// mesh size [cm]
psize = 1.0;

Merge "prolate_mesh_nocap.geo";

Mesh.ColorCarousel = 1;
Mesh.ElementOrder = 1;
Mesh.Lloyd = 1;
//Mesh.HighOrderOptimize = 2;
Mesh.Optimize = 1;
