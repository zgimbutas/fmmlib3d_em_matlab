%
%  Test Maxwell particle FMMs in R^3
%

randn('seed', 1);

nsource = 20
source = randn(3,nsource);

ntarget = nsource;
target = source;
target(3,:) = target(3,:) + 10;

zk = 1.2 + .1*1i;

ifcjvec = 1;
ifcmvec = 1;

cjvec = randn(3,nsource) + 1i*randn(3,nsource);
cmvec = randn(3,nsource) + 1i*randn(3,nsource);

%cjvec = randn(3,nsource);
%cmvec = randn(3,nsource);


ifevec = 1;
ifhvec = 1;
ifevectarg = 1;
ifhvectarg = 1;

'Maxwell dipole, target direct'
[U]=em3dpartdirect_matlab(zk,nsource,source,...
      ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec,...
      ntarget,target,ifevectarg,ifhvectarg);

'Maxwell dipole, target FMM'
iprec = 4;
[V]=emfmm3dpart_matlab(iprec,zk,nsource,source,...
      ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec,...
      ntarget,target,ifevectarg,ifhvectarg);


'Error in H field'
norm(U.hvec - V.hvec)

'Error in E field'
norm(U.evec - V.evec)


'Error in target H field'
norm(U.hvectarg - V.hvectarg)

'Error in target E field'
norm(U.evectarg - V.evectarg)
