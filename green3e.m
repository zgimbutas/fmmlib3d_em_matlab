function evec=green3e(rk,xyz,cjvec);
%GREEN3E Evaluate the dyadic electric Green's function.
%
%   evec=green3e(rk,xyz,cjvec);
%
%   Evaluate the dyadic electric Green's function
%   at the location xyz due to the monochromatic electric dipole
%   cjvec located at the origin.
%
%       ... dyadic Green's function for the electric field 
%       (free space, Lorenz gauge)
%
%                 [        \grad  \grad  ]   exp(I*rk*R)
%       green3e = [ I_3 +  ------------  ]   -----------
%                 [            rk^2      ]        R
%
%       where I_3 is a 3x3 identity matrix
%
%
%   Input parameters:
%
%       rk (complex *16)  - the frequency parameter
%       xyz (real *8) - the target point in R^3
%       cjvec (complex *16) - the strength of the electric dipole   
%
%   Output parameters:
%
%       evec (complex*16) - the value of Green's function at the target
%

n = size(xyz,2);

r = sqrt(sum(xyz.^2,1));

ima = 1i;

fexp = exp(ima.*rk.*r);
f = fexp./r;

%evec = repmat(f,3,1) .* repmat(cjvec,1,n);
evec = kron(f,cjvec);

dx=xyz(1,:);
dy=xyz(2,:);
dz=xyz(3,:);

qmat11 = 2*dx.^2-dy.^2-dz.^2;
qmat22 = 2*dy.^2-dz.^2-dx.^2;
qmat33 = 2*dz.^2-dx.^2-dy.^2;

qmat11 = qmat11.*(1-ima*rk.*r)+(-rk^2.*dx.^2.*r.^2);
qmat22 = qmat22.*(1-ima*rk.*r)+(-rk^2.*dy.^2.*r.^2);
qmat33 = qmat33.*(1-ima*rk.*r)+(-rk^2.*dz.^2.*r.^2);

qmat12 = dx.*dy.*(3-rk^2.*r.^2-3*ima*rk.*r);
qmat23 = dy.*dz.*(3-rk^2.*r.^2-3*ima*rk.*r);
qmat31 = dz.*dx.*(3-rk^2.*r.^2-3*ima*rk.*r);

qmat21 = qmat12;
qmat32 = qmat23;
qmat13 = qmat31;

f = fexp./r.^5/rk^2;

qmat = [qmat11*cjvec(1)+qmat12*cjvec(2)+qmat13*cjvec(3);...
        qmat21*cjvec(1)+qmat22*cjvec(2)+qmat23*cjvec(3);...
        qmat31*cjvec(1)+qmat32*cjvec(2)+qmat33*cjvec(3)];

evec = evec + repmat(f,3,1) .* qmat;

