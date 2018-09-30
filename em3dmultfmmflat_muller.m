function y0=em3dmultfmmflat_muller(A,x0)
%
%     The augmented Muller solver for a dielectric interface
%     with transmission boundary conditions, on a flat triangulated surface. 
%
%     Integral representation is 
%
%     E = (ik) S_k J - grad S_k rho_e + green3m M 
%     H = green3m J - [(ik) S_k M - grad S_k rho_m]
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     Exterior:
%
%        nx (ik) S_k J - nx grad S_k rho_e +(4*pi/2) M + nx grad x S_k M  = nxE
%  +(4*pi/2) rho_e + n. (ik) S_k J - n. grad S_k rho_e + n. grad x S_k M  = n.E
%  +(4*pi/2) J + nx grad x S_k J - [nx (ik) S_k M - nx grad S_k rho_m]    = nxH
% n. grad x S_k J - [+(4*pi/2) rho_m + n. (ik) S_k M - n. grad S_k rho_m] = n.H
%
%     Interior:
%
%        nx (ik) S_k J - nx grad S_k rho_e -(4*pi/2) M + nx grad x S_k M  = nxE
%  -(4*pi/2) rho_e + n. (ik) S_k J - n. grad S_k rho_e + n. grad x S_k M  = n.E
%  -(4*pi/2) J + nx grad x S_k J - [nx (ik) S_k M - nx grad S_k rho_m]    = nxH
% n. grad x S_k J - [-(4*pi/2) rho_m + n. (ik) S_k M - n. grad S_k rho_m] = n.H
%
%     Boundary conditions are
%
%     [nxH] = 0, [nxE] = 0, [n.(eps E)] = 0, [n.(cmu H)] = 0
%
%
disp('em3dmultfmmflat_muller')

x=x0.';

%
% reconstruct the electric and magnetic dipoles from their projections
% on tangential vectors
%
x = reshape(x,6,A.ntri);
jvec = repmat(x(1,:),3,1).*A.triatang1+repmat(x(2,:),3,1).*A.triatang2;
rhoe = x(3,:);
mvec = repmat(x(4,:),3,1).*A.triatang1+repmat(x(5,:),3,1).*A.triatang2;
rhom = x(6,:);

% check, if \int rho_e = 0
rhoe_int=sum(rhoe.*A.triaarea)
% check, if \int rho_m = 0
rhom_int=sum(rhom.*A.triaarea)


%
% evaluate E and H fields on both sides of dielectric interface
%

jvec_e = jvec*A.eps_e;
mvec_e = mvec*A.cmu_e;

[evecj,hvecj]=emfmm3dtria_ccj(A.iprec,A.omega,A.eps_e,A.cmu_e,A.ntri,...
       A.triangles,A.trianorm,A.source,jvec_e,rhoe);

[evecm,hvecm]=emfmm3dtria_ccm(A.iprec,A.omega,A.eps_e,A.cmu_e,A.ntri,...
       A.triangles,A.trianorm,A.source,mvec_e,rhom);

evec_e = evecj + evecm;
hvec_e = hvecj + hvecm;


jvec_i = jvec*A.eps_i;
mvec_i = mvec*A.cmu_i;

[evecj,hvecj]=emfmm3dtria_ccj(A.iprec,A.omega,A.eps_i,A.cmu_i,A.ntri,...
       A.triangles,A.trianorm,A.source,jvec_i,rhoe);

[evecm,hvecm]=emfmm3dtria_ccm(A.iprec,A.omega,A.eps_i,A.cmu_i,A.ntri,...
       A.triangles,A.trianorm,A.source,mvec_i,rhom);

evec_i = evecj + evecm;
hvec_i = hvecj + hvecm;



y = zeros(6,A.ntri);


% [nxH]

c = cross(A.trianorm,hvec_e-hvec_i);

% apply the diagonal term

% interior MFIE
%%%c = c - (4*pi)/2 * (A.eps_i + A.eps_e) * jvec;
% exterior MFIE
c = c + (4*pi)/2 * (A.eps_i + A.eps_e) * jvec;

% project on tangential vectors

y(1,:) = dot(A.triatang1,c);
y(2,:) = dot(A.triatang2,c);


% [nxE]

c = cross(A.trianorm,evec_e-evec_i);

% apply the diagonal term

% interior MFIE
%%%c = c - (4*pi)/2 * (A.cmu_i + A.cmu_e) * mvec;
% exterior MFIE
c = c + (4*pi)/2 * (A.cmu_i + A.cmu_e) * mvec;

% project on tangential vectors

y(4,:) = dot(A.triatang1,c);
y(5,:) = dot(A.triatang2,c);


% [n.(eps E)]

c = dot(A.trianorm,A.eps_e*evec_e-A.eps_i*evec_i);

% apply the diagonal term

% interior EFIE
%%%c = c - (4*pi)/2 * (A.eps_i + A.eps_e) * rhoe;
% exterior EFIE
c = c + (4*pi)/2 * (A.eps_i + A.eps_e) * rhoe;

y(3,:) = c;


% [n.(cmu H)]

c = dot(A.trianorm,A.cmu_e*hvec_e-A.cmu_i*hvec_i);

% apply the diagonal term

% interior EFIE
%%%c = c + (4*pi)/2 * (A.cmu_i + A.cmu_e) * rhom;
% exterior EFIE
c = c - (4*pi)/2 * (A.cmu_i + A.cmu_e) * rhom;

y(6,:) = c;

% flip sign to accelerate the GMRES solver
y(6,:) = -y(6,:);


% optional
% mean zero condition: \int rho_e = 0 (follows from the continuity condition)

%%%y(3,:) = y(3,:) + rhoe_int*ones(1,A.ntri);

% optional
% mean zero condition: \int rho_m = 0 (follows from the continuity condition)

%%%y(6,:) = y(6,:) + rhom_int*ones(1,A.ntri);



y = reshape(y,1,6*A.ntri);

y0=y.';

toc
