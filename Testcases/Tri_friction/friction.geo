Merge "friction.txt";
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(100) = {1};
Mesh.Algorithm = 1;
//+
// Attractors field on points 1 and 55, and on line 1
/*Field[1] = Attractor;*/
/*Field[1].NNodesByEdge = 100; // #attractors on the edges*/
/*Field[1].EdgesList = {616};*/

/*Field[2] = Threshold;*/
/*Field[2].IField = 1;*/
/*Field[2].LcMin = 2;*/
/*Field[2].LcMax = 50;*/
/*Field[2].DistMin = 4;*/
/*Field[2].DistMax = 500;*/

/*Background Field = 50;*/

/*Field[4] = Threshold;*/
/*Field[4].IField = 3;*/
/*Field[4].LcMin = 14;*/
/*Field[4].LcMax = 60;*/
/*Field[4].DistMin = 15;*/
/*Field[4].DistMax = 200;*/
/*//+*/
/*Field[7] = Min;*/
/*Field[7].FieldsList = {2, 4};*/
/*Background Field = 7;*/

/*Mesh.MeshSizeExtendFromBoundary = 0;*/
/*//+*/
/*Mesh.MeshSizeFromCurvature = 0;*/

/*//+*/
/*Print "mesh.png";*/
