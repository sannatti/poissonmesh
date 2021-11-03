function Me = getEdgeInnerProduct(X,Y)

%% Edge inner product in 2D mesh 1/4*sum(P'N^-T V N^-1 P)
%
%  -- 1 -- 
% |       |
% 3       4
% |       |
%  -- 2 -- 
%


m = size(X)-1;
pm = prod(m);

nex = (m(1)+1)*m(2);
ney = (m(2)+1)*m(1);
ne= nex+ney; %number of edge variables

%matrices with order of jx and jy on edges
nexMatrix = reshape(1:nex,m(1)+1,m(2));

neyMatrix = reshape(nex+(1:ney),m(1),m(2)+1);

V = getFaceAreas(X,Y);
[XE1, YE1, XE2, YE2, nEx1, nEy1, nEx2, nEy2] = getEdgeCoord(X,Y);

i=1:m(1); j=1:m(2);

% 1 3
N13 = inv2X2BlockDiagonal(nEx1(i,j),nEy1(i,j),nEx2(i,j),nEy2(i,j));
Ax = nexMatrix(1:m(1),:);
Ay = neyMatrix(:,1:m(2)); 
P13 = sparse(1:2*pm, [Ax(:);Ay(:)],ones(2*pm,1),2*pm,ne);

% 1 4
N14 = inv2X2BlockDiagonal(nEx1(i,j),nEy1(i,j),nEx2(i,j+1),nEy2(i,j+1));
Ax = nexMatrix(1:m(1),:);
Ay = neyMatrix(:,2:m(2)+1);
P14 = sparse(1:2*pm, [Ax(:);Ay(:)],ones(2*pm,1),2*pm,ne);

% 2 3
N23 = inv2X2BlockDiagonal(nEx1(i+1,j),nEy1(i+1,j),nEx2(i,j),nEy2(i,j));
Ax = nexMatrix(2:end,:);
Ay = neyMatrix(:,1:m(2));
P23 = sparse(1:2*pm, [Ax(:);Ay(:)],ones(2*pm,1),2*pm,ne);


% 2 4 
N24 = inv2X2BlockDiagonal(nEx1(i+1,j),nEy1(i+1,j),nEx2(i,j+1),nEy2(i,j+1));
Ax = nexMatrix(2:end,:);
Ay = neyMatrix(:,2:end);
P24 = sparse(1:2*pm, [Ax(:);Ay(:)],ones(2*pm,1),2*pm,ne); 

% 1/4*sum(P'N^-T V N^-1 P)
Me = 0.25* (P13'*N13'*sdiag([V(:);V(:)])*N13*P13 + P14'*N14'*sdiag([V(:);V(:)])*N14*P14+ ...
        P23'*N23'*sdiag([V(:);V(:)])*N23*P23 + P24'*N24'*sdiag([V(:);V(:)])*N24*P24);
    
return
