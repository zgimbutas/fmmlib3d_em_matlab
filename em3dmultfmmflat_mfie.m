function y0=em3dmultfmmflat_mfie(A,x0);
%
%     The exterior MFIE solver with PEC boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     E = green3e J, H = green3m J 
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     -(4*pi/2) J + nx grad x G_k J = nxH  (exterior)
%
%     Boundary conditions:  nxH = J
%
'em3dmultfmmflat_mfie'

x=x0.';

%
% reconstruct the electric dipole from its projection on tangential vectors
%
x = reshape(x,2,A.ntri);
jvec = repmat(x(1,:),3,1).*A.triatang1+repmat(x(2,:),3,1).*A.triatang2;


%
% MFIE operator, H = green3m J = grad x G_k J
%
ifcharge = 1;
ifdipole = 0;
ifpot = 0;
iffld = 1;

dipstr = zeros(1,A.ntri)+1i*zeros(1,A.ntri);

charge = jvec(1,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
y1 = -U.fld;

charge = jvec(2,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
y2 = -U.fld;

charge = jvec(3,:);
[U]=hfmm3dtria(A.iprec,A.zk,A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,charge,ifdipole,dipstr,A.trianorm,ifpot,iffld);
y3 = -U.fld;

% apply the curl operator

hvec=zeros(3,A.ntri)+1i*zeros(3,A.ntri);
hvec(1,:) = y3(2,:)-y2(3,:);
hvec(2,:) = y1(3,:)-y3(1,:);
hvec(3,:) = y2(1,:)-y1(2,:);


% nxH

c = cross(A.trianorm,hvec);

% apply the diagonal term

% interior MFIE
%%%c = c + (4*pi)/2 * jvec;
% exterior MFIE
c = c - (4*pi)/2 * jvec;

% project on tangential vectors

y = zeros(2,A.ntri);
y(1,:) = dot(A.triatang1,c);
y(2,:) = dot(A.triatang2,c);

y = reshape(y,1,2*A.ntri);

y0=y.';

toc
