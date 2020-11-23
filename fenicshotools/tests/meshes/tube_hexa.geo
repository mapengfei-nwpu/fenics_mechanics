Rin  = 0.5;
Rout = 2.0;
Height = 2.0;

Nradius = 4;
Ntheta  = 4;
Nheight = 4;

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

out[] = Extrude
    { 0.0, 0.0, 1.0 }
    { Surface{ surf[] }; Layers{Nheight}; Recombine; };

Physical Volume("Myocardium") = { out[1], out[7], out[13], out[19] };
Physical Surface("Epicardium") = { out[3], out[9], out[15], out[21] };
Physical Surface("Endocardium") = { out[5], out[11], out[17], out[23] };
Physical Surface("Top") = { out[0], out[6], out[12], out[18] };
Physical Surface("Base") = { surf[] };
