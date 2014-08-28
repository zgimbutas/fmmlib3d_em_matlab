%
%  Test Maxwell particle FMMs in R^3
%

randn('seed', 1);

nsource = 20
source = randn(3,nsource);

zk = 1.2 + .1*1i;

ifcjvec = 1;
ifcmvec = 1;

cjvec = randn(3,nsource) + 1i*randn(3,nsource);
cmvec = randn(3,nsource) + 1i*randn(3,nsource);

%cjvec = randn(3,nsource);
%cmvec = randn(3,nsource);


ifevec = 1;
ifhvec = 1;
ifevectarg = 0;
ifhvectarg = 0;

'Maxwell dipole, direct'
[U]=em3dpartdirect_matlab(zk,nsource,source,...
      ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec);

'Maxwell dipole, FMM'
iprec = 4;
[V]=emfmm3dpart_matlab(iprec,zk,nsource,source,...
      ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec);


'Error in H field'
norm(U.hvec - V.hvec)

'Error in E field'
norm(U.evec - V.evec)
