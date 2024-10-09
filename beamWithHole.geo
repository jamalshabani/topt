Point(1) = {0, 0, 0, 0.005};
Point(2) = {1, 0, 0, 0.005};
Point(3) = {1, 1/18, 0, 0.005};
Point(4) = {1, 1/3, 0, 0.005};
Point(5) = {0, 1/3, 0, 0.005};

Point(6) = {4/10, 1/6, 0, 0.005};
Point(7) = {5/10, 1/6, 0, 0.005};
Point(8) = {4/10, 16/60, 0, 0.005};
Point(9) = {3/10, 1/6, 0, 0.005};
Point(10) = {4/10, 4/60, 0, 0.005};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Circle(6) = {10, 6, 7};
Circle(7) = {7, 6, 8};
Circle(8) = {8, 6, 9};
Circle(9) = {9, 6, 10};

Curve Loop(10) = {1, 2, 3, 4, 5};
Curve Loop(11) = {6, 7, 8, 9};

Plane Surface(1) = {10, 11};

Physical Curve("clamped", 7) = {5};
Physical Curve("load", 8) = {2};
Physical Surface("beam", 3) = {1};