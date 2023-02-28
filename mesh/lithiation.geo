// Gmsh project created on Fri Feb 24 15:14:33 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, .0001, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("2", 3) = {1};
//+
Physical Curve("1", 4) = {2};
//+
Physical Surface("5", 5) = {1};
