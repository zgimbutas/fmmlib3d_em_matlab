function hvec=green3m(rk,xyz,cjvec);
%GREEN3M Evaluate the dyadic magnetic Green's function.
%
%   hvec=green3m(rk,xyz,cjvec);
%
%   Evaluate the dyadic magnetic Green's function
%   at the location xyz due to the monochromatic electric dipole
%   cjvec located at the origin.
%
%       ... dyadic Green's function for the magnetic field 
%       (free space, Lorenz gauge)
%
%       green3m =  curl(green3e) = \grad \times green3e
%
%
%                 [        exp(I*rk*R) ]
%       green3m = [ \grad  ----------- ] \times
%                 [             R      ]
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
%       hvec (complex*16) - the value of Green's function at the target
%

n = size(xyz,2);

r = sqrt(sum(xyz.^2,1));

ima = 1i;

fexp = exp(ima.*rk.*r);
f = (1-ima*rk.*r).*fexp./r.^3;

fder = repmat(f,3,1) .* xyz;
hvec = -cross(fder,repmat(cjvec,1,n));
