Rin  = 0.500000;
Rout = 2.000000;
psize = 1.50000;

Point(1) = { Rin,  0.0, 0.0, psize };
Point(2) = { Rout, 0.0, 0.0, psize };

Line(1) = { 1, 2 };

line_id = 1;
For num In {0:2}
    out[] = Extrude
        { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
        { Line{line_id}; };
    line_id = out[0];
    surf[num]    = out[1];
    bc_epi[num]  = out[2];
    bc_endo[num] = out[3];
EndFor
out[] = Extrude
    { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
    { Line{line_id}; };
surf[num]    = out[1];
bc_epi[num]  = out[2];
bc_endo[num] = out[3];

Physical Surface("Myocardium") = { surf[] };
Physical Line("Epicardium")    = { bc_epi[] };
Physical Line("Endocardium")   = { bc_endo[] };

