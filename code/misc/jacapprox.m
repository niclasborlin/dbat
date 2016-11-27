function J=jacapprox(fun,x,h,params)
%JACAPPROX Numerical approximation of jacobian.
%
%   J=JACAPPROX(F,X,H), where F is a function name or handle, X is an
%   N-vector, and H is a scalar, computes a numerical approximation of the
%   Jacobian of the function F at X using a central difference of length 2*H.
%   J is returned as a sparse array. If H is a vector, uses 2*H(I) along
%   dimension I.
%
%   J=JACAPPROX(F,X) uses H=1e-6.
%
%   J=JACAPPROX(F,X,H,P), where P is a cell array, sends P{:} as extra
%   parameters to the function F.


if nargin<3, h=1e-6; end
if nargin<4, params={}; end

% Evaluate F(X).
f0=feval(fun,x,params{:});

% Pre-allocate sparse Jacobian.
n=length(x);
J=sparse(length(f0),n);

if isscalar(h)
    h=repmat(h,n,1);
end

I=speye(n);

for i=1:n
    % Numerical approximation with central difference
    %      F( X+H(I)*I(:,I) ) - F( X-H(I)*I(:,I) ) 
    %      ---------------------------------------
    %                       2*H(I)
    J(:,i)=(feval(fun,x+h(i)*I(:,i),params{:})-...
            feval(fun,x-h(i)*I(:,i),params{:}))/(2*h(i));
end
