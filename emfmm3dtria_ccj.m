function [evec,hvec]=emfmm3dtria_ccj(iprec,omega,eps,cmu,ntri,triaflat,trianorm,source,jvec,rhoe)
%
%     Integral representation is 
%
%     E = (i omega cmu) G_k J - grad G_k rho_e, H = green3m J 
%

zk=omega*sqrt(eps)*sqrt(cmu);

%
% MFIE operator, H = green3m J = grad x G_k J
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
iffld = 1;

dipstr = zeros(1,ntri)+1i*zeros(1,ntri);

charge = jvec(1,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x1 = U.pot;
y1 = -U.fld;

charge = jvec(2,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x2 = U.pot;
y2 = -U.fld;

charge = jvec(3,:);
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);
x3 = U.pot;
y3 = -U.fld;

% apply the curl operator

hvec=zeros(3,ntri)+1i*zeros(3,ntri);
hvec(1,:) = y3(2,:)-y2(3,:);
hvec(2,:) = y1(3,:)-y3(1,:);
hvec(3,:) = y2(1,:)-y1(2,:);



%
% EFIE operator, E = (i omega cmu) G_k J - grad G_k rho_e
%
evec=zeros(3,ntri)+1i*zeros(3,ntri);
evec(1,:) = x1;
evec(2,:) = x2;
evec(3,:) = x3;

evec = evec*1i*omega*cmu;

ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 1;

charge = rhoe;
[U]=hfmm3dtria(iprec,zk,ntri,triaflat,trianorm,source,...
       ifcharge,charge,ifdipole,dipstr,trianorm,ifpot,iffld);

evec = evec + U.fld;


