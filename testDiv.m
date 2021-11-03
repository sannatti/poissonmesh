%% Derivative test for Div operator
% equation used for testing: u(x,y)= exp(-sig*exp(x^2 + y^2)



for k = 2
n = 2^k;

%% generate the grid mesh 
[X,Y] = ndgrid(linspace(-1,1,n), linspace(-1,1,n));
m = size(X)-1;

% nodal points of the mesh
c = 0.2;
X = X + c*sin(pi*X).*sin(pi*Y);
Y = Y + c*sin(pi*X).*sin(pi*Y);


% 
[XC,YC]                                         = getCellCoord(X,Y);
[XE1, YE1, XE2, YE2, nEx1, nEy1, nEx2, nEy2]    = getEdgeCoord(X,Y); 
L                                               = getEdgeLenghts(X,Y);
area                                            = getFaceAreas(X,Y);

%get DIV operator
DIVinit = getDivergence(X,Y);
DIV = sdiag(1./area(:))*DIVinit*sdiag(L);

%function test u
sig = 1;
u   = @(X,Y) exp(-sig*(X.^2 + Y.^2) );

%the first derivative functions
ux  = @(X,Y) -sig*2*X.*u(X,Y);
uy  = @(X,Y) -sig*2*Y.*u(X,Y);


%the second derivative functions
uxx = @(X,Y) -sig*2*(u(X,Y) + X.*ux(X,Y));
uxy = @(X,Y) -sig*2*X.*uy(X,Y);
uyx = @(X,Y) uxy(X,Y);
uyy = @(X,Y) -sig*2*(u(X,Y) + Y.*uy(X,Y));

% gradient of u  
JF1 = ux(XE1(:),YE1(:)).*nEx1(:) + ...
      uy(XE1(:),YE1(:)).*nEy1(:);
JF2 = ux(XE2(:),YE2(:)).*nEx2(:) + ...
      uy(XE2(:),YE2(:)).*nEy2(:);



%the know divergence in analytical form 
q   = uxx(XC,YC) +  uyy(XC,YC);

% the test function on cell coordinates
uu  = u(XC,YC);

% the gradient on edges
J   = [JF1(:); JF2(:)];


% numerical and analytical difference DIV
r1 = DIV*J - q(:);

fprintf('error = %3.2e, mesh size n = %3.2e \n', max(abs(r1)),n)
%%%%%%

%test function
Fx = @(X,Y) 2*cos(pi*X) + X.*Y.^2;
Fy = @(X,Y) 1+X*0 + 3*cos(pi*Y);

%partial derivatives
Fxx = @(X,Y) -pi*2.*sin(pi*X)+ Y.^2 +X*0;
Fyy = @(X,Y) -pi*3.*sin(pi*Y);

Fx1 = Fx(XE1, YE1);
Fy1 = Fy(XE1,YE1);
Fx2 = Fx(XE2,YE2);
Fy2 = Fy(XE2, YE2);



F = [Fx1(:).*nEx1(:)+Fy1(:).*nEy1(:);Fx2(:).*nEx2(:)+Fy2(:).*nEy2(:)];
Q = Fxx(XC,YC) + Fyy(XC,YC);

err = DIV*F-Q(:);
 
 
%fprintf('2 error = %3.2e, n = %3.2e \n', max(abs(err)),n)
end 


figure(1)
subplot(3,1,1)
surf(reshape((err),m)); colorbar
title('error')
subplot(3,1,2)
surf(reshape(Q,m)); colorbar
title('analytical function')
subplot(3,1,3)
surf(reshape(DIV*F,m)); colorbar
title('numerical function') 

figure(2)
subplot(3,1,1)
surf(reshape(abs(r1),m)); colorbar
title('|error|')
subplot(3,1,2)
surf(reshape(abs(q),m)); colorbar
title('|analytical|')
subplot(3,1,3)
surf(reshape(abs(DIV*J),m)); colorbar
title('|numerical|')


figure(3)
plot(DIV*J,'r')
hold on
plot(q(:),'b')
legend('numerical','analytical')
title('numerical and analytical divergence')
hold off





