// Gmsh project created on Mon Dec 20 05:18:46 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.2};
//+
Point(2) = {0, 1, 0, 0.2};
//+
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
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
