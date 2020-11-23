psize = 0.5;
base = 1.0;
height = 3.0;

// Center
Point(1) = { 0.0, 0.0, 0.0 };

// Endocardium
Point(2) = { base, 0.0, 0.0, psize };
Point(3) = { 0.0, height, 0.0, psize };
Ellipse(1) = { 3, 1, 3, 2 };

Extrude
    { { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, Pi/3}
    { Line{1}; }

