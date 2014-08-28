%
%  Test Maxwell triangle FMMs in R^3
%

%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere180.a.tri';
filename_geo = 'sphere320.a.tri';
%filename_geo = 'sphere1280.a.tri';
%filename_geo = 'sphere2880.a.tri';
%filename_geo = 'sphere11520.a.tri';
%filename_geo = 'sphere20480.a.tri';
%filename_geo = 'cube_r1.a.tri';
%filename_geo = 'cube_r2.a.tri';
%filename_geo = 'cube_r3.a.tri';
%filename_geo = 'cube_r4.a.tri';
%filename_geo = '../hellskitchen/Data/cube768.a.tri';
%filename_geo = '../hellskitchen/Data/cube3072.a.tri';
%filename_geo = '../hellskitchen/Data/cube12288.a.tri';

[verts,ifaces,nverts,nfaces] = atriread(filename_geo);
nverts,nfaces

%  adjust refined cube
%verts = verts-0.5

%  ellipsoid
%verts(3,:) = verts(3,:)*1.5;

%
%  Construct triangle vertex, normal, area, and centroid arrays
%

ntri = nfaces;
[triangles,trianorm,triaarea,source,triatang1,triatang2]=...
     atriproc(verts,ifaces);

%plot3(source(1,:),source(2,:),source(3,:),'*')

nsource = ntri;

%
%  Initialize matrix multiplication routine
%

% EFIE and MFIE resonances

%zk= .8733491486*pi;
%zk= 1.231935202*pi;
%zk= 1.583088865*pi;
%zk= 1.929578411*pi;
%zk= 2.272804959*pi;
%zk= 2.613592175*pi;
%zk= 2.952471726*pi;

%zk= 1.430296653*pi;
%zk= 1.834566041*pi;
%zk= 2.224327839*pi;
%zk= 2.604590204*pi;
%zk= 2.978047487*pi;


%%zk=2*pi/6
zk=1
A.zk = zk;

nsource = ntri;


A.iprec = 0;

A.ntri = ntri;
A.triangles = triangles;
A.trianorm = trianorm;
A.triatang1 = triatang1;
A.triatang2 = triatang2;
A.triaarea = triaarea;
A.source = source;


%
%  Construct the test right hand side
%
%%%rhs = ones(3,nsource);

ntest_source = 1;
test_source = [0.01;-0.03;0.02];
test_source = [100;-200;300];

iftest_cjvec = 1;
iftest_cmvec = 0;
%%test_cjvec = [1;-3;2];
%%test_cmvec = [1;-2;3];
test_cjvec = [1;1;1];
test_cmvec = [1;1;1];

ifevec = 0;
ifhvec = 0;
ifevectarg = 1;
ifhvectarg = 1;

[U]=em3dpartdirecttarg_matlab(zk,ntest_source,test_source,...
      iftest_cjvec,test_cjvec,iftest_cmvec,test_cmvec,...
      ifevec,ifhvec,nsource,source,ifevectarg,ifhvectarg);

%%%rhs = U.hvectarg;

nxe = cross(A.trianorm,U.evectarg);
rhs = zeros(3,A.ntri);
rhs(1,:) = dot(A.triatang1,nxe);
rhs(2,:) = dot(A.triatang2,nxe);
rhs(3,:) = dot(A.trianorm,U.evectarg);

M = sqrt(triaarea)';
M = repmat(M,1,3);
M = reshape(M,1,3*nsource)';
%
%  Call the solver
%
'Interior/Exterior augmented EFIE solver for the Maxwell equations in R^3'
tic
rhs0 = reshape(rhs,1,3*nsource);
sol0 = gmres_simple(@(x) em3dmultfmmflat_aefie(A,x), rhs0.', 1e-2, 20);
%%%sol0 = gmres_simple(@(x) M.*em3dmultfmmflat_aefie(A,x./M), M.*rhs0.', 1e-2, 20)./M;
sol0 = sol0.';
sol = reshape(sol0,3,nsource);
time_gmres=toc

jvec = repmat(sol(1,:),3,1).*A.triatang1+repmat(sol(2,:),3,1).*A.triatang2;
rhoe = sol(3,:);

%
%  Evaluate fields directly
%
ntarget = 1;
target = [100;-200;300];
target = [0.01;-0.03;0.02];

%%%target = [100;200;-300];

ifevec = 0;
ifhvec = 0;
ifevectarg = 1;
ifhvectarg = 1;
[U]=em3dpartdirecttarg_matlab(zk,ntest_source,test_source,...
      iftest_cjvec,test_cjvec,iftest_cmvec,test_cmvec,...
      ifevec,ifhvec,ntarget,target,ifevectarg,ifhvectarg);

if( ifevectarg == 1 ), evec = U.evectarg; end
if( ifhvectarg == 1 ), hvec = U.hvectarg; end


%
%  Evaluate fields via the solution of integral equation
%
tic

ifcjvec=1;
ifcmvec=0;
cjvec=jvec.*repmat(triaarea,3,1);
rho_e=rhoe.*repmat(triaarea,1,1);
cmvec=zeros(3,nsource);
rho_m=zeros(1,nsource);

if( 2 == 2 ),
[U]=em3dpartdirecttarg_matlab(zk,nsource,source,...
      ifcjvec,cjvec,ifcmvec,cmvec,...
      ifevec,ifhvec,ntarget,target,ifevectarg,ifhvectarg);
else
[U]=em3dccpartdirecttarg(zk,nsource,source,...
      ifcjvec,cjvec,rho_e,ifcmvec,cmvec,rho_m,...
      ifevec,ifhvec,ntarget,target,ifevectarg,ifhvectarg);
end

toc

%
%  Finally, print the results
%
'Fields at the targets'

if( ifevectarg == 1 ),
evec,U.evectarg
rel_error_evec=((evec-U.evectarg)./evec)
end

if( ifhvectarg == 1 ),
hvec,U.hvectarg
rel_error_hvec=((hvec-U.hvectarg)./hvec)
end
