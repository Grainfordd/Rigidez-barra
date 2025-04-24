Geometry.PointLabels = 1;
tm = 1e6;
Point(1) = {0, 0, 0, tm};
Point(2) = {0, 10, 0, tm};
Point(3) = {10, 0, 0, tm};
Point(4) = {10, 10, 0, tm};

Line(1) = {1, 3};
Line(2) = {1, 4};
Line(3) = {1, 2};
Line(4) = {2, 4};
Line(5) = {2, 3};
Line(6) = {3, 4};

// Physical Line("Estructura")  = {1,2, 3, 4, 5, 6};
Physical Line("Estructura") = { Line{:} };
Physical Point("Fijo") = {4};
Physical Point("Rod_x") = {3};


Mesh.MshFileVersion = 2.2;
Mesh 1;
Save "malla.msh";
