function [U]=emfmm3dpart_matlab(iprec,zk,nsource,source,ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec,ntarget,target,ifevectarg,ifhvectarg)
%EMFMM3DPART Maxwell dipole FMM in R^3.
%
% Maxwell FMM in R^3: evaluate all pairwise dipole
% interactions (ignoring self-interaction) and interactions with targets.
%
%
% [U]=EMFMM3DPART(IPREC,ZK,NSOURCE,SOURCE,...
%         IFCJVEC,CJVEC,IFCMVEC,CMVEC,IFEVEC,IFHVEC);
%
% [U]=EMFMM3DPART(IPREC,ZK,NSOURCE,SOURCE,...
%         IFCJVEC,CJVEC,IFCMVEC,CMVEC,IFEVEC,IFHVEC,...
%         NTARGET,TARGET);
%
% [U]=EMFMM3DPART(IPREC,ZK,NSOURCE,SOURCE,...
%         IFCJVEC,CJVEC,IFCMVEC,CMVEC,IFEVEC,IFHVEC,...
%         NTARGET,TARGET,IFEVECTARG,IFHVECTARG);
%
%
% This subroutine evaluates the Maxwell potential and field due
% to a collection of electric and magnetic dipoles. 
% We use (exp(ikr)/r) for the Green's function, without the (1/4 pi) scaling. 
% Self-interactions are not-included.
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% zk - complex Helmholtz parameter
% nsource - number of sources
% source - real (3,nsource): source locations
% ifcjvec - electric dipole computation flag
%
%         0 => do not compute
%         1 => include electric dipole contribution
% 
% cjvec - complex (3,nsource): electric dipole strengths 
%
% ifehvec - (E,H) computation flag, 
%         1 => compute the electric and magnetic fields, otherwise no
%
% ntarget - number of targets
% target - real (3,ntarget): target locations
%
% ifehvectarg - target (E,H) computation flag, 
%         1 => compute the target electric and magnetic fields, otherwise no
%
%
% Output parameters: 
%
% U.evec - complex (3,nsource) - E field at source locations
% U.hvec - complex (3,nsource) - H field at source locations
%
% U.evectarg - complex (3,target) - E field at target locations
% U.hvectarg - complex (3,target) - H field at target locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%

if( ifevec ~= ifhvec ),
  error('emfmm3dpart: not implemented')
end

if( ifevectarg ~= ifhvectarg ),
  error('emfmm3dpart: not implemented')
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifevectarg = 0;
  ifhvectarg = 0;
end

if( nargin == 12 ) 
  ifevectarg = 0;
  ifhvectarg = 0;
end


ifcjvec = double(ifcjvec); ifcmvec = double(ifcmvec); 

ifevec = double(ifevec); ifhvec = double(ifhvec); 
ifevectarg = double(ifevec); ifhvectarg = double(ifhvectarg); 

evec=zeros(3,1);
hvec=zeros(3,1);
evectarg=zeros(3,1);
hvectarg=zeros(3,1);

if( ifevec == 1 ), evec=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifhvec == 1 ), hvec=zeros(3,nsource)+1i*zeros(3,nsource); end;

if( ifevectarg == 1 ), evectarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;
if( ifhvectarg == 1 ), hvectarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;


ier=0;

if( ifcjvec == 1 ),
%
% MFIE operator, H = green3m J = grad x G_k J
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
iffld = 1;

dipstr = zeros(1,nsource)+1i*zeros(1,nsource);
dipvec = zeros(3,nsource);

charge = cjvec(1,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x1 = U.pot;
y1 = -U.fld;

charge = cjvec(2,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x2 = U.pot;
y2 = -U.fld;

charge = cjvec(3,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x3 = U.pot;
y3 = -U.fld;

% apply the curl operator

hvec(1,:) = hvec(1,:) + (y3(2,:)-y2(3,:));
hvec(2,:) = hvec(2,:) + (y1(3,:)-y3(1,:));
hvec(3,:) = hvec(3,:) + (y2(1,:)-y1(2,:));


%
% EFIE operator, E = (i k) green3e J
% EFIE operator, E = (i k) G_k J + 1/ (i k )  grad grad G_k J 
%
evec(1,:) = evec(1,:) + x1*1i*zk;
evec(2,:) = evec(2,:) + x2*1i*zk;
evec(3,:) = evec(3,:) + x3*1i*zk;

ifcharge = 0;
ifdipole = 1;
ifpot = 0;
iffld = 1;

charge = zeros(1,nsource);

dipstr = cjvec(1,:);
dipvec = [ones(1,nsource); zeros(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

evec = evec - U.fld / (1i*zk);

dipstr = cjvec(2,:);
dipvec = [zeros(1,nsource); ones(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

evec = evec - U.fld / (1i*zk);

dipstr = cjvec(3,:);
dipvec = [zeros(1,nsource); zeros(1,nsource); ones(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

evec = evec - U.fld / (1i*zk);


end


if( ifcmvec == 1 ),
%
% MFIE operator, E = green3m M = grad x G_k M
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
iffld = 1;

dipstr = zeros(1,nsource)+1i*zeros(1,nsource);
dipvec = zeros(3,nsource);

charge = cmvec(1,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x1 = U.pot;
y1 = -U.fld;

charge = cmvec(2,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x2 = U.pot;
y2 = -U.fld;

charge = cmvec(3,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
x3 = U.pot;
y3 = -U.fld;

% apply the curl operator

evec(1,:) = evec(1,:) + (y3(2,:)-y2(3,:));
evec(2,:) = evec(2,:) + (y1(3,:)-y3(1,:));
evec(3,:) = evec(3,:) + (y2(1,:)-y1(2,:));


%
% EFIE operator, H = -(i k) green3e M
% EFIE operator, H = -(i k) G_k M - 1/ (i k )  grad grad G_k M
%
hvec(1,:) = hvec(1,:) - x1*1i*zk;
hvec(2,:) = hvec(2,:) - x2*1i*zk;
hvec(3,:) = hvec(3,:) - x3*1i*zk;

ifcharge = 0;
ifdipole = 1;
ifpot = 0;
iffld = 1;

dipstr = cmvec(1,:);
dipvec = [ones(1,nsource); zeros(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

hvec = hvec + U.fld / (1i*zk);

dipstr = cmvec(2,:);
dipvec = [zeros(1,nsource); ones(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

hvec = hvec + U.fld / (1i*zk);

dipstr = cmvec(3,:);
dipvec = [zeros(1,nsource); zeros(1,nsource); ones(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);

hvec = hvec + U.fld / (1i*zk);

end



if( ifcjvec == 1 && ntarget > 0 ),
%
% MFIE operator, H = green3m J = grad x G_k J
%
ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 0;
ifpottarg = 1;
iffldtarg = 1;

dipstr = zeros(1,nsource)+1i*zeros(1,nsource);
dipvec = zeros(3,nsource);

charge = cjvec(1,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);
x1 = U.pottarg;
y1 = -U.fldtarg;

charge = cjvec(2,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);
x2 = U.pottarg;
y2 = -U.fldtarg;

charge = cjvec(3,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);
x3 = U.pottarg;
y3 = -U.fldtarg;

% apply the curl operator

hvectarg(1,:) = hvectarg(1,:) + (y3(2,:)-y2(3,:));
hvectarg(2,:) = hvectarg(2,:) + (y1(3,:)-y3(1,:));
hvectarg(3,:) = hvectarg(3,:) + (y2(1,:)-y1(2,:));


%
% EFIE operator, E = (i k) green3e J
% EFIE operator, E = (i k) G_k J + 1/ (i k )  grad grad G_k J 
%
evectarg(1,:) = evectarg(1,:) + x1*1i*zk;
evectarg(2,:) = evectarg(2,:) + x2*1i*zk;
evectarg(3,:) = evectarg(3,:) + x3*1i*zk;

ifcharge = 0;
ifdipole = 1;
ifpot = 0;
iffld = 0;
ifpottarg = 0;
iffldtarg = 1;

charge = zeros(1,nsource);

dipstr = cjvec(1,:);
dipvec = [ones(1,nsource); zeros(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);

evectarg = evectarg - U.fldtarg / (1i*zk);

dipstr = cjvec(2,:);
dipvec = [zeros(1,nsource); ones(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);

evectarg = evectarg - U.fldtarg / (1i*zk);

dipstr = cjvec(3,:);
dipvec = [zeros(1,nsource); zeros(1,nsource); ones(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,
       ntarget,target,ifpottarg,iffldtarg);

evectarg = evectarg - U.fldtarg / (1i*zk);


end


if( ifcmvec == 1 && ntarget > 0 ),
%
% MFIE operator, E = green3m M = grad x G_k M
%
ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 0;
ifpottarg = 1;
iffldtarg = 1;

dipstr = zeros(1,nsource)+1i*zeros(1,nsource);
dipvec = zeros(3,nsource);

charge = cmvec(1,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);
x1 = U.pottarg;
y1 = -U.fldtarg;

charge = cmvec(2,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);
x2 = U.pottarg;
y2 = -U.fldtarg;

charge = cmvec(3,:);
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);
x3 = U.pottarg;
y3 = -U.fldtarg;

% apply the curl operator

evectarg(1,:) = evectarg(1,:) + (y3(2,:)-y2(3,:));
evectarg(2,:) = evectarg(2,:) + (y1(3,:)-y3(1,:));
evectarg(3,:) = evectarg(3,:) + (y2(1,:)-y1(2,:));


%
% EFIE operator, H = -(i k) green3e M
% EFIE operator, H = -(i k) G_k M - 1/ (i k )  grad grad G_k M
%
hvectarg(1,:) = hvectarg(1,:) - x1*1i*zk;
hvectarg(2,:) = hvectarg(2,:) - x2*1i*zk;
hvectarg(3,:) = hvectarg(3,:) - x3*1i*zk;

ifcharge = 0;
ifdipole = 1;
ifpot = 0;
iffld = 0;
ifpottarg = 0;
iffldtarg = 1;

dipstr = cmvec(1,:);
dipvec = [ones(1,nsource); zeros(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);

hvectarg = hvectarg + U.fldtarg / (1i*zk);

dipstr = cmvec(2,:);
dipvec = [zeros(1,nsource); ones(1,nsource); zeros(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);

hvectarg = hvectarg + U.fldtarg / (1i*zk);

dipstr = cmvec(3,:);
dipvec = [zeros(1,nsource); zeros(1,nsource); ones(1,nsource)];
[U]=hfmm3dpart(iprec,zk,nsource,source,...
       ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,...
       ntarget,target,ifpottarg,iffldtarg);

hvectarg = hvectarg + U.fldtarg / (1i*zk);

end






if( ifevec == 1 ), U.evec=evec; end
if( ifhvec == 1 ), U.hvec=hvec; end
if( ifevectarg == 1 ), U.evectarg=evectarg; end
if( ifhvectarg == 1 ), U.hvectarg=hvectarg; end
U.ier=ier;


