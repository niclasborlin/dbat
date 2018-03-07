function J=jacapprox(fun,x,h,returnFull)
%JACAPPROX Numerical approximation of jacobian.
%
%   J=JACAPPROX(F,X,H), where F is a function name or handle, X is an
%   N-vector, and H is a scalar, computes a numerical approximation of the
%   Jacobian of the function F at X using a central difference of length 2*H.
%   J is returned as a sparse array. If H is a vector, uses 2*H(I) along
%   dimension I.
% 
%   If X is an array, it is unrolled before calling F. If F returns
%   an array, it is unrolled before being used in the computation.
%
%   J=JACAPPROX(F,X) uses H=1e-6.
%
%   J is typically returned as a sparse matrix. Use
%   J=JACAPPROX(F,X,H,TRUE) to force the return value to be full.

switch nargin
  case {0,1}
    narginchk(2,4);
  case 2
    h=1e-6;
    returnFull=false;
  case 3
    if islogical(h)
        returnFull=h;
        h=1e-6;
    else
        returnFull=false;
    end
  case 4
    % Do nothing
end

x=x(:);

% Evaluate F(X).
f0=vec(feval(fun,x));

% Pre-allocate Jacobian.
n=length(x);
if returnFull
    J=nan(length(f0),n);
    I=eye(n);
else
    J=sparse(length(f0),n);
    I=speye(n);
end

if isscalar(h)
    h=repmat(h,n,1);
end

for i=1:n
    % Numerical approximation with central difference
    %      F( X+H(I)*I(:,I) ) - F( X-H(I)*I(:,I) ) 
    %      ---------------------------------------
    %                       2*H(I)
    J(:,i)=vec((feval(fun,x+h(i)*I(:,i))-...
                feval(fun,x-h(i)*I(:,i)))/(2*h(i)));
end


function v=vec(x)
% Specify vec here to make jacapprox self-contained.

v=x(:);
