Point(1) = {0,  0, 0, 0.005};
Point(2) = {1, 0, 0, 0.005};
Point(3) = {1, 1/6, 0, 0.005};
Point(4) = {1,  1, 0, 0.005};
Point(5) = {0,  1, 0, 0.005};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Curve Loop(6) = {1, 2, 3, 4, 5};

Plane Surface(1) = {6};

Physical Curve("left", 7) = {5};
Physical Curve("rightcorner", 8) = {2};
Physical Surface("beam", 2) = {1};
