function [A,dA,dAn]=affine2mat(b,varargin)
%AFFINE2MAT 2D affine transformation matrix.
%
%   A=AFFINE2MAT(B) returns the affine transformation matrix
%   [1+B(1), B(2); 0, 1].
%
%   [A,dA]=... also returns a struct dA with the analytical Jacobian
%   with respect to B in the field dB. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(b), A=selftest(nargin>1 && varargin{1}); return; end

% Otherwise, verify number of parameters.
narginchk(1,2);

A=[]; %#ok<NASGU>
dA=[];
dAn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dA=struct('dB',[]);
    dAn=dA;
end

% What Jacobians to compute?
cB=nargout>1 && (length(varargin)<1 || varargin{1});

%% Test parameters
if any(size(b)~=[2,1])
    error([mfilename,': bad size']);
end

%% Actual function code
A=[1+b(1),b(2);0,1];

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cB
        fun=@(b)feval(mfilename,b);
        dAn.dB=jacapprox(fun,b);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cB
        dA.dB=sparse([1,3],[1,2],1,4,2);
    end
end


function fail=selftest(verbose)

% Set up test data.
b=rand(2,1)+1;

fail=full_self_test(mfilename,{b},1e-8,1e-8,verbose);
