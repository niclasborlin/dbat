function [P,PP,res]=pm_resect_3pt(X,x,use,behind,relax)
%PM_RESECT_3PT Spatial resection from three points.
%
%[P,PP,res]=pm_resect_3pt(X,x,use[,behind[,relax]])
%X      - 3xn or 4xn matrix of 3D object points in world coordinate system.
%x      - 2xn or 3xn matrix of 2D measured points in normalized camera
%         coordinates.
%use    - 3-index vector (or boolean) of which points to use for
%         resection. Remaining points are used to distinguish between
%         possible 3-pt solutions.
%behind - are the points behind the camera? Default: false.
%relax  - include complex polynomial solutions. Default: false.
%P      - best normalized camera matrix, or [] if none can be found.
%PP     - m-cell array of all suggested camera matrices.
%res    - m-vector with average projection residuals for each camera.
%
%   References:
%     Haralick, Lee, Ottenberg, Nölle (1994), "Review and Analysis of Solutions
%       of the Three Point Perspective Pose Estimation Problem".
%       Internation Journal of Computer Vision, 13(3):331-356.
%     McGlone, Mikhail, Bethel, eds. (2004), "Manual of Photogrammetry",
%       5th ed., Chapter 11.1.3.4, pp. 786-788. American Society of
%       Photogrammetry and Remote Sensing.


if nargin<4, behind=false; end
if nargin<5, relax=false; end

% Convert points to euclidean.
if size(X,1)==4, X=euclidean(X); end
if size(x,1)==3, x=euclidean(x); end

if nnz(use)~=3
    error('Can only use 3 points for resection');
end

% Extract points to use for calculation.
XTest=X;
xTest=x;
X=X(:,use);
x=x(:,use);

% Append focal length and normalize to get direction vector corresponding
% to each projected point.
x=[x;ones(1,size(x,2))];
x=x./repmat(sqrt(sum(x.^2)),3,1);

% 3D distances.
a=norm(diff(X(:,[2,3]),[],2));
b=norm(diff(X(:,[1,3]),[],2));
c=norm(diff(X(:,[1,2]),[],2));

% Angles between direction vectors.
alpha=subspace(x(:,2),x(:,3));
beta=subspace(x(:,1),x(:,3));
gamma=subspace(x(:,1),x(:,2));

% Set up 4th order polynomial.
a2mc2db2=(a^2-c^2)/b^2;
a2pc2db2=(a^2+c^2)/b^2;
b2mc2db2=(b^2-c^2)/b^2;
b2ma2db2=(b^2-a^2)/b^2;

A4=(a2mc2db2-1)^2-4*c^2/b^2*cos(alpha)^2;
A3=4*(a2mc2db2*(1-a2mc2db2)*cos(beta)+2*c^2/b^2*cos(alpha)^2*cos(beta)-(1-a2pc2db2)*cos(alpha)*cos(gamma));
A2=2*(a2mc2db2^2+2*a2mc2db2^2*cos(beta)^2+2*b2mc2db2*cos(alpha)^2+2*b2ma2db2*cos(gamma)^2-4*a2pc2db2*cos(alpha)*cos(beta)*cos(gamma)-1);
A1=4*(-a2mc2db2*(1+a2mc2db2)*cos(beta)+2*a^2/b^2*cos(gamma)^2*cos(beta)-(1-a2pc2db2)*cos(alpha)*cos(gamma));
A0=(1+a2mc2db2)^2-4*a^2/b^2*cos(gamma)^2;

% Extract roots.
v=roots([A4,A3,A2,A1,A0]);
if ~relax
    % Allow only real roots.
    v=v(abs(imag(v)./abs(v))<1e-3);
    v=real(v);
else
    % Accept all solutions. Testing will determine which one is best.
    v=unique(real(v));
end

% Corresponding u values.
u=((-1+a2mc2db2)*v.^2-2*a2mc2db2*cos(beta)*v+1+a2mc2db2)./...
  (2*(cos(gamma)-v*cos(alpha)));

% Calculate distances from u, v.
s12=b^2./(1+v.^2-2*v*cos(beta));
s1=sqrt(s12);
s3=v.*s1;
s2=u.*s1;

% Only keep unique positive distances.
valid=s1>=0 & s2>=0 & s3>=0;
s1=s1(valid);
s2=s2(valid);
s3=s3(valid);
s123=unique([s1,s2,s3],'rows');

% Calculate rigid body transformation for each solution.
res=zeros(1,size(s123,1));
PP=cell(size(res));
for i=1:size(s123,1)
    % Absolute orientation (McGlone et al., Ch. 11.1.3.4.2).
    
    % 3D points in camera coordinate system.
    cxi=repmat(s123(i,:),3,1).*x;
    
    if behind
        cxi=-cxi;
    end
    
    ob=X(:,3)-X(:,1);
    oc=X(:,2)-X(:,1);
    cb=cxi(:,3)-cxi(:,1);
    cc=cxi(:,2)-cxi(:,1);

    r1=ob/norm(ob);
    r2=cross(ob,oc);
    r2=r2/norm(r2);
    r3=cross(ob,cross(ob,oc));
    r3=r3/norm(r3);
    
    oRdelta=[r1,r2,r3];

    r1=cb/norm(cb);
    r2=cross(cb,cc);
    r2=r2/norm(r2);
    r3=cross(cb,cross(cb,cc));
    r3=r3/norm(r3);
    cRdelta=[r1,r2,r3];

    cRo=cRdelta*oRdelta';

    oxO=X(:,1)-cRo'*cxi(:,1);
    
    P=cRo*[eye(3),-oxO];
    PP{i}=P;
    % Calculate mean pt residual for projected test points.
    res(i)=sqrt(mean(sum((euclidean(P*homogeneous(XTest))-xTest).^2)));
end

[~,i]=min(res);
if (~isempty(i))
    P=PP{i(1)};
else
    P=[];
end

