function [A,dA,dAn]=affine2mat(b1,b2,varargin)
%AFFINE2MAT 2D affine transformation matrix.
%
%   A=AFFINE2MAT(B1,B2) returns the affine transformation matrix
%   [1+B1, B2; 0, 1].
%
%   [A,dA]=... also returns a struct dA with the analytical Jacobian
%   with respect to B in the field dB. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(b1), A=selftest(nargin>1 && b2); return; end

% Otherwise, verify number of parameters.
narginchk(2,4);

A=[];
dA=[];
dAn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dA=struct('dB1',[],...
              'dB2',[]);
    dAn=dA;
end

% What Jacobians to compute?
cB1=nargout>1 && (length(varargin)<1 || varargin{1});
cB2=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
if ~isscalar(b1) || ~isscalar(b2)
    error([mfilename,': bad size']);
end

%% Actual function code
A=[1+b1,b2;0,1];

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cB1
        fun=@(b1)feval(mfilename,b1,b2);
        dAn.dB1=jacapprox(fun,b1);
    end
    if cB2
        fun=@(b2)feval(mfilename,b1,b2);
        dAn.dB2=jacapprox(fun,b2);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cB1
        dA.dB1=sparse(1,1,1,4,1);
    end
    if cB2
        dA.dB2=sparse(3,1,1,4,1);
    end
end


function fail=selftest(verbose)

% Set up test data.
b1=1+rand;
b2=1+rand;

fail=full_self_test(mfilename,{b1,b2},1e-8,1e-8,verbose);
