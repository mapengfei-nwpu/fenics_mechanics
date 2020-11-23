Rin  = 0.500000;
Rout = 2.000000;
Nradius = 4;
Ntheta  = 4;


Point(1) = { Rin,  0.0, 0.0 };
Point(2) = { Rout, 0.0, 0.0 };

Line(1) = { 1, 2 };
Transfinite Line{1} = Nradius;

line_id = 1;
For num In {0:2}
    out[] = Extrude
        { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
        { Line{line_id}; Layers{Ntheta}; Recombine; };
    line_id = out[0];
    surf[num]    = out[1];
    bc_epi[num]  = out[2];
    bc_endo[num] = out[3];
EndFor
out[] = Extrude
    { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
    { Line{line_id}; Layers{Ntheta}; Recombine; };
surf[num]    = out[1];
bc_epi[num]  = out[2];
bc_endo[num] = out[3];

Physical Surface("Myocardium") = { surf[] };
Physical Line("Epicardium")    = { bc_epi[] };
Physical Line("Endocardium")   = { bc_endo[] };
