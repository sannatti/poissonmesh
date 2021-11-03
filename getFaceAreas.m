function area = getFaceAreas(X,Y)

%% Face areas
m = size(X)-1;

% A----B 
% |    |
% D----C
%
% A = |AC x BD|*0.5

% compute the diagonals

i=1:m(1); j=1:m(2);
dxp = X(i+1,j+1)-X(i,j);
dyp = Y(i+1,j+1)-Y(i,j);

dxm = X(i,j+1)-X(i+1,j);
dym = Y(i,j+1)-Y(i+1,j);

% the cross product 
crossp = dxp.*dym - dyp.*dxm;
normn1 = sqrt(crossp.^2);
area   = normn1/2;