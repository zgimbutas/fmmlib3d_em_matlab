function y0=em3dmultfmmflat_aefie(A,x0);
%
%     The exterior augmented EFIE solver with PEC boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     E = (ik) G_k J - grad G_k rho_e, H = green3m J 
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%                       nx (ik) G_k J - nx grad G_k rho_e = nxE 
%     -(4*pi/2) rho_e + n. (ik) G_k J - n. grad G_k rho_e = n.E (exterior)
%
%     Boundary conditions:  nxE = 0, n.E = rho_e
%
'em3dmultfmmflat_aefie'

x=x0.';

%
% reconstruct the electric dipole from its projection on tangential vectors
%
x = reshape(x,3,A.ntri);
jvec = repmat(x(1,:),3,1).*A.triatang1+repmat(x(2,:),3,1).*A.triatang2;
rhoe = x(3,:);

% check, if \int rho_e = 0
rhoe_int=sum(rhoe.*A.triaarea)

%
% EFIE operator, E = (ik) G_k J - grad G_k rho_e
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
iffld = 0;

dipstr = zeros(1,A.ntri)+1i*zeros(1,A.ntri);

charge = jvec(1,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
x1 = U.pot;

charge = jvec(2,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
x2 = U.pot;

charge = jvec(3,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
x3 = U.pot;


evec=zeros(3,A.ntri)+1i*zeros(3,A.ntri);
evec(1,:) = x1;
evec(2,:) = x2;
evec(3,:) = x3;

evec = evec*1i*A.zk;


ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 1;

charge = rhoe;
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);

evec = evec + U.fld;


% nxE

c = cross(A.trianorm,evec);

% project on tangential vectors

y = zeros(3,A.ntri);
y(1,:) = dot(A.triatang1,c);
y(2,:) = dot(A.triatang2,c);

% n.E

c = dot(A.trianorm,evec);

% apply the diagonal term

% interior EFIE
%%%c = c + (4*pi)/2 * rhoe;
% exterior EFIE
c = c - (4*pi)/2 * rhoe;

y(3,:) = c;


y = reshape(y,1,3*A.ntri);

y0=y.';

toc
