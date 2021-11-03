function [XE1, YE1, XE2, YE2, nEx1, nEy1, nEx2, nEy2] = getEdgeCoord(X,Y)
% Edge coordinates of 2D mesh
% Input: 
% X, Y nodal coordainates
% Output: 

%number of cells
m = size(X)-1; 



%% Edge coordinates
%first direction (down)
i = 1:m(1)+1; j = 1:m(2);
XE1 = 1/2*(X(i,j+1)+ X(i,j)); 
YE1 = 1/2*(Y(i,j+1)+ Y(i,j));
%second direction (horizontal)
i = 1:m(1); j = 1:m(2)+1;
XE2 = 1/2*(X(i,j)+ X(i+1,j));
YE2 = 1/2*(Y(i,j)+ Y(i+1,j));


%% Edge tangents
i = 1:m(1)+1; j = 1:m(2);
tEx1 = X(i,j+1)-X(i,j); 
tEy1 = Y(i,j+1)-Y(i,j);
% normalize edge tangents
norm1 = sqrt(tEx1.^2+tEy1.^2);
tEx1 = tEx1./norm1;
tEy1 = tEy1./norm1;

i = 1:m(1); j = 1:m(2)+1;
tEx2 = X(i+1,j)-X(i,j); 
tEy2 = Y(i+1,j)-Y(i,j);
% normalize edge tangents
norm2 = sqrt(tEx2.^2+tEy2.^2);
tEx2 = tEx2./norm2;
tEy2 = tEy2./norm2;


%% Edge normals
nEx2 = -tEy2;
nEy2 = tEx2; 

nEx1 = tEy1;
nEy1 = -tEx1;
