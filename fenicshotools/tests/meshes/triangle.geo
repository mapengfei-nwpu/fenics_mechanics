psize = 2.0;

Point(1) = { 0.0, 0.0, 0.0, psize };
Point(2) = { 1.0, 0.0, 0.0, psize };
Point(3) = { 0.0, 1.0, 0.0, psize };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 1 };

Line Loop(1) = { 3, 1, 2 };
Plane Surface(1) = { 1 };

