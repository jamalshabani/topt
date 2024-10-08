Point(1) = {0, 0, 0, 0.005};
Point(2) = {1, 0, 0, 0.005};
Point(3) = {1, 1/3, 0, 0.005};
Point(4) = {0, 1/3, 0, 0.005};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(5) = {1, 2, 3, 4};

Plane Surface(1) = {5};

Physical Curve("clamped", 7) = {1};
Physical Curve("load", 8) = {3};
Physical Surface("beam", 2) = {1};