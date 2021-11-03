function DIV = getDivergence(X,Y)
%% DIV = getDivergence(X,Y)
% 2D Div: from edges to cell centers

dx = @(n)(spdiags(ones(n+1,1)*[-1,1],[0,1],n,n+1));

m = size(X)-1;

DIV = [kron(speye(m(2)),dx(m(1))) , ...
      kron(dx(m(2)),speye(m(1)))];
  
 

