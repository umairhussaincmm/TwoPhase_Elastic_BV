// Gmsh project created on Tue Feb 28 10:04:10 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
//Point{1} In Surface {1};
//+
Physical Point("1", 2) = {1};
//+
Physical Curve("2", 3) = {1};
//+
Physical Surface("5", 4) = {1};
