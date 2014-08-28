function [evec,hvec]=em3dipole3mt(rk,source,target,cmvec)
%EM3DIPOLE3MT: Evaluate E and H fields of a magnetic dipole.
%
%   [EVEC,HVEC]=em3dipole3mt(RK,SOURCE,TARGET,CMVEC);
%
%   Evaluate E and H fields at the location target due
%   to the monochromatic magnetic dipole cmvec located
%   at an arbitrary source location.
%
%   Input parameters:
%
%       rk (complex *16)  - the frequency parameter
%       source (real *8 ) - the source point in R^3
%       target (real *8 ) - the target point in R^3
%       cmvec (complex *16) - the strength of the magnetic dipole
%
%   Output parameters:
%
%       evec (complex*16) - the electric field at the target
%       hvec (complex*16) - the magnetic field at the target
%

n = size(target,2);
xyz = target - repmat(source,1,n);

evec = green3m(rk,xyz,cmvec);
hvec = green3e(rk,xyz,cmvec);

ima = 1i;
hvec = hvec * (-ima*rk);
