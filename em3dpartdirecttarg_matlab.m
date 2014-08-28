function [U]=em3dpartdirecttarg_matlab(zk,nsource,source,ifcjvec,cjvec,ifcmvec,cmvec,ifevec,ifhvec,ntarget,target,ifevectarg,ifhvectarg)
%EM3DPARTDIRECTTARG Maxwell dipole interactions in R^3, direct evaluation.
%
% Maxwell interactions in R^3: evaluate all pairwise dipole
% interactions (ignoring self-interaction) and interactions with targets.
%


if( nargin == 7 ) 
  ifevec = 1;
  ifhvec = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifevectarg = 0;
  ifhvectarg = 0;
end

if( nargin == 9 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifevectarg = 0;
  ifhvectarg = 0;
end

if( nargin == 11 ) 
  ifevectarg = 1;
  ifhvectarg = 1;
end

evec=0;
hvec=0;
evectarg=0;
hvectarg=0;

if( ifevec == 1 ), evec=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifhvec == 1 ), hvec=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifevectarg == 1 ), evectarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;
if( ifhvectarg == 1 ), hvectarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;

for j=1:nsource
for i=1:nsource
  if( i ~= j ), 
    if( ifcjvec ),
    [evec0,hvec0] = em3dipole3et(zk,source(:,i),source(:,j),cjvec(:,i));
    if( ifevec ), evec(:,j) = evec(:,j) + evec0; end
    if( ifhvec ), hvec(:,j) = hvec(:,j) + hvec0; end
    end
    if( ifcmvec ),
    [evec0,hvec0] = em3dipole3mt(zk,source(:,i),source(:,j),cmvec(:,i));
    if( ifevec ), evec(:,j) = evec(:,j) + evec0; end
    if( ifhvec ), hvec(:,j) = hvec(:,j) + hvec0; end
    end
  end
end
end

for j=1:ntarget
for i=1:nsource
  if( ifcjvec ),
  [evectarg0,hvectarg0] = em3dipole3et(zk,source(:,i),target(:,j),cjvec(:,i));
  if( ifevectarg ), evectarg(:,j) = evectarg(:,j) + evectarg0; end
  if( ifhvectarg ), hvectarg(:,j) = hvectarg(:,j) + hvectarg0; end
  end
  if( ifcmvec ),
  [evectarg0,hvectarg0] = em3dipole3mt(zk,source(:,i),target(:,j),cmvec(:,i));
  if( ifevectarg ), evectarg(:,j) = evectarg(:,j) + evectarg0; end
  if( ifhvectarg ), hvectarg(:,j) = hvectarg(:,j) + hvectarg0; end
  end
end
end


if( ifevec == 1 ), U.evec=evec; end;
if( ifhvec == 1 ), U.hvec=hvec; end;
if( ifevectarg == 1 ), U.evectarg=evectarg; end;
if( ifhvectarg == 1 ), U.hvectarg=hvectarg; end;



