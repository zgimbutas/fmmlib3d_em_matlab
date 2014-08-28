function [evec,hvec]=em3dipole3et(rk,source,target,cjvec)
%EM3DIPOLE3ET: Evaluate E and H fields of an electric dipole.
%
%   [EVEC,HVEC]=em3dipole3et(RK,SOURCE,TARGET,CJVEC);
%
%   Evaluate E and H fields at the location target due
%   to the monochromatic electric dipole cjvec located 
%   at an arbitrary source location.
%
%   Input parameters:
%
%       rk (complex *16)  - the frequency parameter
%       source (real *8 ) - the source point in R^3
%       target (real *8 ) - the target point in R^3
%       cjvec (complex *16) - the strength of the electric dipole   
%
%   Output parameters:
%
%       evec (complex*16) - the electric field at the target
%       hvec (complex*16) - the magnetic field at the target
%

n = size(target,2);
xyz = target - repmat(source,1,n);

evec = green3e(rk,xyz,cjvec);
hvec = green3m(rk,xyz,cjvec);

ima = 1i;
evec = evec * (ima*rk);
