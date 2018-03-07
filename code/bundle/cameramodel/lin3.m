function [Q,dQ,dQn]=lin3(P,c,varargin)
%LIN3 3D scaling for the DBAT projection model
%
%   Q=LIN3(P,C) lins the 3D points in the 3-by-N array P by the
%   scalar C.
%
%   [Q,dQ]=... also returns a struct dQ with the analytical Jacobians
%   with respect to P and C in the fields dP and dC. For more
%   details, see DBAT_BUNDLE_FUNCTIONS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

if nargin>=1 && ischar(P), Q=selftest(nargin>1 && c); return; end

Q=[];
dQ=[];
dQn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dQ=struct('dP',[],...
              'dC',[]);
    dQn=dQ;
end

% What Jacobians to compute?
cP=nargout>1 && (length(varargin)<1 || varargin{1});
cC=nargout>1 && (length(varargin)<2 || varargin{2});

%% Test parameters
[m,n]=size(P);
if m~=3 || ~isscalar(c)
    error([mfilename,': bad size']);
end

%% Actual function code
Q=c*P;

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    if cP
        fmt=@(P)reshape(P,3,[]);
        fun=@(P)feval(mfilename,fmt(P),c);
        dQn.dP=jacapprox(fun,P);
    end
    if cC
        fun=@(c)feval(mfilename,P,c);
        dQn.dC=jacapprox(fun,c);
    end
end

if nargout>1
    %% Analytical Jacobian
    if cP
        dQ.dP=c*speye(numel(P));
    end
    if cC
        dQ.dC=vec(P);
    end    
end


function fail=selftest(verbose)

% Set up test data.
n=3;
m=5;
c=rand+1;
P=rand(n,m);

% Compute Jacobians.
[~,A,N]=feval(mfilename,P,c);

fail=comparejacobianstructs(mfilename,A,N,1e-8,1e-8,verbose);
