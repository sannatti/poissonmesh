function L = getEdgeLenghts(X,Y)

m = size(X)-1;

%% Edge lenghts
Lx1 = Y(:,2:end) - Y(:,1:m(2));
Lx2 = X(:,2:end) - X(:,1:m(2));
Lx  = sqrt(Lx1.^2 + Lx2.^2);
Ly1 = X(2:end,:) - X(1:m(1),:);
Ly2 = Y(2:end,:) - Y(1:m(1),:);
Ly  = sqrt(Ly1.^2 +Ly2.^2);
L   = [Lx(:);Ly(:)];
L   = abs(L);


