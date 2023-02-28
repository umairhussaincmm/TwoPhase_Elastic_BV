// Gmsh project created on Wed Apr 21 05:37:49 2021
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 1.0};
//+
Point(3) = {5, 3, 0, 1.0};
//+
Point(4) = {0, 3, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("1") = {4, 2};
//+
Physical Curve("2") = {3, 1};
//+
Physical Surface("5") = {1};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Curve {4, 2} = 30 Using Progression 1;
//+
Transfinite Curve {1, 3} = 50 Using Progression 1;
