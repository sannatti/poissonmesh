function A = sdiag(a)

A = spdiags(a(:),0,numel(a),numel(a));