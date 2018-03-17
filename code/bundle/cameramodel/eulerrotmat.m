function [M,dM,dMn]=eulerrotmat(ang,seq,fixed)
%EULERROTMAT 3D Euler rotation matrix.
%
%   M=EULERROTMAT(A,SEQ,FIXED) computes the 3D rotation matrix that
%   correspond to the Euler angles in the 3-vector A for the axis
%   sequence given by the 3-digit integer SEQ. Each digit can be 1
%   (X), 2 (Y), or 3 (Z), and determines the axes about which the
%   elementary rotations takes place. For instance, SEQ=123
%   corresponds to x-y-z rotations, SEQ=313 corresponds to z-x-z
%   rotations. If FIXED is TRUE, the rotations are performed in a
%   fixed coordinate system. If FIXED is FALSE, the rotations are
%   performed in a moving coordinate system.
%
%   [M,dM]=... also returns a struct dM with the analytical Jacobian
%   with respect to K in the field dA. For more details, see
%   DBAT_BUNDLE_FUNCTIONS.
%
%   EULERROTMAT is defined for any angles and any sequence of axis.
%   However, if two subsequent axis are identical, or if the second
%   rotation rotates the first axis to the third, the decomposition
%   from rotation matrix to angles is not unique.
%
%   References:
%
%     Lucas (1963), "Differentiation of the Orientation Matrix by Matrix
%       Multipliers". Photogrammetric Engineering 29(4):708-715.
%
%     Forstner, Wrobel (2004), "Mathematical Concepts in Photogrammetry",  
%       Ch. 2.1.2. In "Manual of Photogrammetry", 5th ed. McGlone, et al.,
%       eds. ASPRS.
%
%SEE ALSO: DBAT_BUNDLE_FUNCTIONS

% Treat selftest call separately.
if nargin>=1 && ischar(ang), M=selftest(nargin>1 && seq); return; end

% Otherwise, verify number of parameters.
narginchk(3,3);

M=[]; %#ok<NASGU>
dM=[];
dMn=[];

if nargout>1
    % Construct empty Jacobian struct.
    dM=struct('dA',[]);
    dMn=dM;
end

%% Test parameters
if length(ang)~=3
    error([mfilename,': bad size']);
end

%% Actual function code

% Rotation functions.
fh={@R1,@R2,@R3};

% Decode what rotations we should use.
i1=floor(seq/100);
i2=floor(rem(seq,100)/10);
i3=rem(seq,10);

% Compute the elementary rotations.
if nargout<2
    M1=fh{i1}(ang(1));
    M2=fh{i2}(ang(2));
    M3=fh{i3}(ang(3));
else
    [M1,P1]=fh{i1}(ang(1));
    [M2,P2]=fh{i2}(ang(2));
    [M3,P3]=fh{i3}(ang(3));
end

% Combine differently depending on whether the rotation is with
% respect to a fixed or rotation frame.
if fixed
    M=M3*M2*M1;
else
    M=M1*M2*M3;
end

if nargout>2
    %% Numerical Jacobian

    % FMT is function handle to repackage vector argument to what
    % the function expects.
    fun=@(ang)feval(mfilename,ang,seq,fixed);
    dMn.dA=jacapprox(fun,ang);
end

if nargout>1
    %% Analytical Jacobian
    if fixed
        dA1=M*P1;
        dA2=M3*M2*P2*M1;
        dA3=P3*M;
    else
        dA1=P1*M;
        dA2=M1*M2*P2*M3;
        dA3=M*P3;
    end
    dM.dA=[dA1(:),dA2(:),dA3(:)];
end

% Elementary rotations about each axis.

function [R,P]=R1(alpha)

R=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
P=[0,0,0;0,0,-1;0,1,0];


function [R,P]=R2(alpha)

R=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];
P=[0,0,1;0,0,0;-1,0,0];


function [R,P]=R3(alpha)

R=[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
P=[0,-1,0;1,0,0;0,0,0];


function fail=selftest(verbose)

% Set up test data.
seqs=[];
for i=1:3
    for j=1:3 % find(~ismember(1:3,i))
        for k=1:3 % find(~ismember(1:3,j))
            seqs(end+1)=i*100+j*10+k; %#ok<AGROW>
        end
    end
end

ang=rand(3,1);
fail=false;
for i=1:length(seqs)
    seq=seqs(i);
    for t=[false,true]
        fail=fail | full_self_test(mfilename,{ang,seq,t},1e-8,1e-8,verbose,0);
    end
end
