function B = inv2X2BlockDiagonal(a11,a12,a21,a22)

% inverse matrix of 2x2

a11=a11(:); a12=a12(:); 
a21=a21(:); a22=a22(:);

detA = a11.*a22-a12.*a21;


b11 =  a22./detA;
b12 = -a12./detA;
b21 = -a21./detA;
b22 = a11./detA;

B = [ 
    sdiag(b11) sdiag(b12) 
    sdiag(b21) sdiag(b22)
    ];


detB = b11.*b22-b12.*b21;
return;