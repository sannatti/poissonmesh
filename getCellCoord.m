function [XC, YC] = getCellCoord(X,Y)
%% [XC, YC] = getCellCoord(X,Y)
% in 2D mesh

m = size(X)-1;

i = 1:m(1); j = 1:m(2);
XC = 1/4*(X(i,j)+ X(i+1,j) + X(i,j+1) + X(i+1,j+1));
YC = 1/4*(Y(i,j)+ Y(i+1,j) + Y(i,j+1) + Y(i+1,j+1));