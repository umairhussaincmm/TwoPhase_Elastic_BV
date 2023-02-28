// Gmsh project created on Tue Jan  3 11:22:42 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.2};
//+
Point(2) = {0, 5e-6, 0, 0.2};
//+
Circle(1) = {0, 0, 0, 5e-6, 0, 2*Pi};
//+
Physical Point("1") = {2};
//+
Physical Curve("2") = {1};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("5") = {1};
