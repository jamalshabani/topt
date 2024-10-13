Point(1) = {0, 0, 0, 0.05};
Point(2) = {1, 0, 0, 0.05};
Point(3) = {0, 1, 0, 0.05};
Point(4) = {1, 1, 0, 0.05};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {2, 4};
Line(5) = {4, 3};

Curve Loop(1) = {1, 2, 3};
Curve Loop(2) = {-2, 4, 5};
Plane Surface(1) = {1};
Plane Surface(2) = {2};