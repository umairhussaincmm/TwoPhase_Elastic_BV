//+
Point(1) = {0, 0, 0, .5};
//+
Point(2) = {70, 0, 0, .5};
//+
Point(3) = {0, 70, 0, .5};
//+
Point(4) = {-70, 0, 0, .5};
//+
Point(5) = {0, -70, 0, .5};
//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 70, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("1") = {1};
//+
Physical Surface("5") = {1};

