function [x,ier]=bicgstab_simple(A,b,tol,maxit,x0)
% BiCG(stab) algorithm
%
%  This subroutine solves a complex linear system Ax=b by means
%  of BiCG(stab) algorithm.
%
%  X = BICGSTAB_SIMPLE(A,b);
%  X = BICGSTAB_SIMPLE(A,b,tol);
%  X = BICGSTAB_SIMPLE(A,b,tol,maxit);
%  X = BICGSTAB_SIMPLE(A,b,tol,maxit,x0);
%
%  Input parameters:
%
%  A - the matrix of the system (or the function handle to evaluate A(x) )
%  b - the right hand side 
%  tol - the required accuracy
%  maxit - the maximum number of iteration permitted
%  x0 - the initial guess
%
%  Output parameters:
%
%  x - the solution of the system
%  ier - error return code
%     ier=0 normal execution of the subroutine
%     ier=4 means that the maximum number iterations maxit
%           has been reached without achieving the required accuracy tol
%


[n,m]=size(b);

if( nargin < 3 ) tol = 1e-6; end
if( nargin < 4 ) maxit = min(n,20); end

ier=0;

% A is a matrix, construct function handle 
if( isnumeric(A) ), A = (@(x) A*x); end

% initialize 
if( nargin == 5 ),
xk=x0;
else
xk=zeros(n,1);
end
axk=A(xk);

if( n == 1 )
  xk = b/axk;
  axk = b;
  return;
end


% find the first direction
er = 0;
er1 = 0;

ek=b-axk;
ekm1 = ek;

bek = ek;
bekm1 = bek;

norm1=sqrt(ek'*ek);
er1=sqrt(bek'*ek);



% bicgstab algorithm
for i=2:maxit

aekm1=A(ekm1);

cd1=bek'*ek;
cd2=bek'*aekm1;

alpha=cd1/cd2;

esk=ek-alpha*aekm1;
etk=A(esk);

cd1=etk'*esk;
cd2=etk'*etk;

omega=cd1/cd2;

% update the solution
xk=xk+alpha*ekm1+omega*esk;

% update the residual
ek=esk-omega*etk;

er=sqrt(ek'*ek);
er2=sqrt(bek'*ek);

beta=(er2/er1)*(alpha/omega);
er1=er2;

ekm1=ek+beta*(ekm1-omega*aekm1);

%%%norm(r,2)
%%%fprintf('iter: %d, norm(r,2): %f\n',i,norm(r,2));
%%%fprintf('iter: %d, rms=norm(r,2)/sqrt(n): %17.12f\n',i,er/sqrt(n));
fprintf('iter: %d, rel=norm(r,2)/norm(b,2): %17.12e\n',i,er/norm(b,2));

%%%if( i == 1 ) er0=er; end;
if( i == 1 ) norm1=norm(b,2); end;
if( er < tol*norm1 ), x=xk; return; end;

end

% % no convergence, abort, return the solution
x=xk;
ier=4;
