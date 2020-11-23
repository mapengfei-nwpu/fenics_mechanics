// ========================
// generates a prolate mesh
// ------------------------

// x = d * sinh(lambda) * sin(mu) * cos(theta)
// y = d * sinh(lambda) * sin(mu) * sin(theta)
// z = d * cosh(lambda) * cos(mu)

// Note:
// - mesh is translated in order to have the base at z = 0
// - mu_epi such that z_endo = z_epi at the base

// Things
mu_base_endo = mu_base / 180.0 * Pi;
mu_base_epi  = Acos(Cosh(l_endo)/Cosh(l_epi) * Cos(mu_base_endo));

mu_apex = 0.0;
theta = 0.0;

x_endo_base = dfocal * Sinh(l_endo) * Sin(mu_base_endo) * Cos(theta);
y_endo_base = dfocal * Cosh(l_endo) * Cos(mu_base_endo) * Sin(theta);
z_endo_base = dfocal * Cosh(l_endo) * Cos(mu_base_endo);
x_endo_apex = dfocal * Sinh(l_endo) * Sin(mu_apex) * Cos(theta);
y_endo_apex = dfocal * Cosh(l_endo) * Cos(mu_apex) * Sin(theta);
z_endo_apex = dfocal * Cosh(l_endo) * Cos(mu_apex);

x_epi_base = dfocal * Sinh(l_epi) * Sin(mu_base_epi) * Cos(theta);
y_epi_base = dfocal * Cosh(l_epi) * Cos(mu_base_epi) * Sin(theta);
z_epi_base = dfocal * Cosh(l_epi) * Cos(mu_base_epi);
x_epi_apex = dfocal * Sinh(l_epi) * Sin(mu_apex) * Cos(theta);
y_epi_apex = dfocal * Cosh(l_epi) * Cos(mu_apex) * Sin(theta);
z_epi_apex = dfocal * Cosh(l_epi) * Cos(mu_apex);

// Center
Point(1) = { 0.0, 0.0, - z_endo_base };

// Endocardium
Point(2) = { x_endo_base, y_endo_base, 0.0, psize };
Point(3) = { x_endo_apex, y_endo_apex, z_endo_apex - z_endo_base, psize };
Ellipse(1) = { 3, 1, 3, 2 };

// Epicardium
Point(4) = { x_epi_base, y_epi_base, z_epi_base - z_endo_base, psize };
Point(5) = { x_epi_apex, y_epi_apex, z_epi_apex - z_endo_base, psize };
Ellipse(2) = { 5, 1, 5, 4 };

// Base
Line(3) = { 2, 4 };

// Apex
Line(4) = { 3, 5 };

// Section
Line Loop(1) = { 2, -3, -1, 4 };
Plane Surface(1) = { 1 };

// Extrusion
surf_id = 1;
For num In {0:3}
  out[] = Extrude
    { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2}
    { Surface{surf_id}; };
  vol[num] = out[1];
  surf_epi[num] = out[2];
  surf_base[num] = out[3];
  surf_endo[num] = out[4];
  surf_id = out[0];
EndFor

// Physical entities
// -----------------
Physical Volume(1) = { vol[] };
// Base
Physical Surface(10) = { surf_base[] };
// Epi
Physical Surface(20) = { surf_epi[] };
// Endo
Physical Surface(30) = { surf_endo[] };

