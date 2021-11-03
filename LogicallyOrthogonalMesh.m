
% Derviative test for numerical solver of Poisson's equation: 
% Div(f) = u, where f = 0 on boudary, 
% on logically orthogonal mesh
% Sanna Tyrvainen 2014


%grid size
for k = 3:7
 n = 2^k; % number of grid points

[X,Y] = ndgrid(linspace(-1,1,n),linspace(-1,1,n/2));
m = size(X)-1;

%create mesh with following nodal points
c = 0.1; % c = 0 produces regular mesh 
X = X + c*sin(pi*X).*sin(pi*Y);
Y = Y + c*sin(pi*X).*sin(pi*Y);

% Edge lenghts
L = getEdgeLenghts(X,Y);
% Face Areas
area = getFaceAreas(X,Y);
% Edge coordinates
[XE1, YE1, XE2, YE2, nEx1, nEy1, nEx2, nEy2] = getEdgeCoord(X,Y);
% Cell center coordinates
[XC,YC] = getCellCoord(X,Y);

%% Get Divergence matrix
DIV = getDivergence(X,Y); 


% Edge inner product 1/4*sum(P'N^-T V N^-1 P)
Me = getEdgeInnerProduct(X,Y);
empty = zeros(prod(m),prod(m));

A = [Me', sdiag(L)*DIV'; -DIV*sdiag(L) empty]; 

    
%% a test problem Div(f) = q with f=0 on boundary

u = @(x,y) sin(pi*x).*sin(pi*y);
J = @(x,y) pi*cos(pi*x).*sin(pi*y);
q = @(x,y) -2*pi^2*u(x,y);

uAnalytic = u(XC,YC);   % Analytic u at cell centres;
Jx1 = J(XE1,YE1);      % grad U x direction
Jy1 = J(YE1,XE1);      % grad U x direction
Jx2 = J(XE2,YE2);
Jy2 = J(YE2,XE2);        % grad U y direction
J = [Jx1(:).*nEx1(:)+ Jy1(:).*nEy1(:); Jx2(:).*nEx2(:)+Jy2(:).*nEy2(:)];

RHSNumeric = A*[J;uAnalytic(:)];
qAnalytic = q(XC,YC);
e = zeros(size(J,1),1);

qNumeric = -RHSNumeric(size(e,1)+1:end);


LHS = [e(:);- sdiag(area(:))*qAnalytic(:)];
ju = A\LHS;
uNumeric = ju(size(e,1)+1:end);
JNumeric = ju(1:size(e,1));

error1 = max(abs(qAnalytic(:)-sdiag(1./area(:))*qNumeric(:)));
error2 = max(abs(uAnalytic(:)-uNumeric(:)));
error3 = max(abs(J(:)-JNumeric(:)));


fprintf('%3.0f &  %3.2e & %3.2e  &  %3.2e\n', n, error2, error1, error3)


end 


figure(1)
subplot(3,1,1)
plot(((J(:)- JNumeric(:)))); 
title('derivative error')
subplot(3,1,2)
plot((J)); 
title('derivative analytic')
subplot(3,1,3)
plot((JNumeric));
title('derivative numeric')



figure(2)
subplot(3,1,1)
surf(reshape((qAnalytic(:)-qNumeric(:)),m)); colorbar
title('q error')
subplot(3,1,2)
surf(reshape(qAnalytic,m)); colorbar
title('q analytic')
subplot(3,1,3)
surf(reshape(qNumeric,m)); colorbar
title('q numeric')

figure(3)
subplot(3,1,1)
surf(reshape((uAnalytic(:)-uNumeric(:)),m)); colorbar
title('u error')
subplot(3,1,2)
surf(reshape(uAnalytic,m)); colorbar
title('u numeric')
subplot(3,1,3)
surf(reshape(uNumeric,m)); colorbar
title('u analytic')


figure(4)
plot(X,Y,'r-x',X',Y','r-x')
hold on
%quiver(XE1,YE1,nEx1,nEy1) % edge normals
%quiver(XE2,YE2,nEx2,nEy2) 
title('the mesh')
hold off

