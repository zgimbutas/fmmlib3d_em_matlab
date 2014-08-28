function [evec,hvec]=emfmm3dtria_ccm(iprec,omega,eps,cmu,ntri,triaflat,trianorm,source,mvec,rhom)
%
%     Integral representation is 
%
%     E = green3m M, H = - [(i omega eps) G_k M - grad G_k rho_m]
%

zk=omega*sqrt(eps)*sqrt(cmu);

%
% dual MFIE operator, E = green3m M = grad x G_k M
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
iffld = 1;

dipstr = zeros(1,ntri)+1i*zeros(1,ntri);

charge = mvec(1,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x1 = U.pot;
y1 = -U.fld;

charge = mvec(2,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x2 = U.pot;
y2 = -U.fld;

charge = mvec(3,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x3 = U.pot;
y3 = -U.fld;

% apply the curl operator

evec=zeros(3,ntri)+1i*zeros(3,ntri);
evec(1,:) = y3(2,:)-y2(3,:);
evec(2,:) = y1(3,:)-y3(1,:);
evec(3,:) = y2(1,:)-y1(2,:);



%
% dual EFIE operator, H = - [(i omega eps) G_k M - grad G_k rho_m]
%
hvec=zeros(3,ntri)+1i*zeros(3,ntri);
hvec(1,:) = x1;
hvec(2,:) = x2;
hvec(3,:) = x3;

hvec = hvec*1i*omega*eps;

ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 1;

charge = rhom;
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);

hvec = hvec + U.fld;

hvec = -hvec;
