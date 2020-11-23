psize = 2.0;

Point(1) = { 0.0, 0.0, 0.0, psize };
Point(2) = { 1.0, 0.0, 0.0, psize };
Point(3) = { 0.0, 1.0, 0.0, psize };
Point(4) = { 0.0, 0.0, 1.0, psize };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 1, 3 };
Line(4) = { 1, 4 };
Line(5) = { 3, 4 };
Line(6) = { 2, 4 };

Line Loop(7) = {6, -4, 1};
Plane Surface(8) = {7};
Line Loop(9) = {2, -3, 1};
Plane Surface(10) = {9};
Line Loop(11) = {5, -4, 3};
Plane Surface(12) = {11};
Line Loop(13) = {2, 5, -6};
Plane Surface(14) = {13};
Surface Loop(15) = {14, 10, 12, 8};
Volume(16) = {15};
