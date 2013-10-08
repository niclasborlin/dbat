function v=pm_eulerbundle1p_v3(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,PI,constr,orthoP,cIO,cEO,cOP,cPI,cSI,varargin)
%PM_EULERBUNDLE1_V Veto function for Euler cam bundle adjustment.
%
%v=pm_eulerbundle1_v3(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%or
%v=pm_eulerbundle1_v3('setup',[dMin,dMax],ltr,upDir,OPix,OPid,camID)
%x      - vector of unknowns.
%pt     - measured projected points.
%ptCams - physical camera for each point.
%IO     - vector with camera inner orientation as columns
%         pp - principal point [xp;yp].
%         f  - focal length.
%         K  - radial lens distortion coefficients.
%         P  - tangential lens distortion coefficients.
%         a  - with affine lens distortion coefficients.
%         u  - image unit.
%nK     - number of radial koefficients.
%nP     - number of tangential koefficients.
%EO     - matrix with photo camera outer orientation as columns
%         C     - camera center [Xc,Yc,Zc].
%         ang   - Euler angles.
%         ax    - Euler axis sequence
%                 0 - 'xyz'.
%                 1 - 'zxz'.
%cams   - vector of camera numbers for each photo.
%OP     - 3 x nObj matrix with object coordinates
%vis    - sparse matrix indicating if obj point i is visible in photo j.
%cIO    - logical matrix indicating which elements of IO are unknown.
%cEO    - logical matrix indicating which elements of EO are unknown.
%cOP    - logical matrix indicating which elements of OP are unknown.
%dMin,
%dMax   - depth limits in project units.
%ltr    - camera indices showing the required left-to-right sequence.
%upDir  - 'up' direction for the definition of left and right.
%OPix   - object point indices to use for calculation center of object.
%OPid   - object point ids.
%camId  - camera ids.
%v      - true if any object points are behind the corresponding cameras
%         or any camera pairs are in the wrong order.

persistent DEPTH_LIMITS UPDIR LTR OPIX OPID VETO_ACTIVATED CAMID

if ischar(x)
    if x(1)=='s'
        % ...('setup',[dMin,dMax],ltr,upDir,OPix,OPid)
        DEPTH_LIMITS=pt;
        LTR=ptCams;
        if ~iscell(LTR)
            LTR={LTR};
        end
        UPDIR=IO;
        OPIX=nK;
        OPID=nP;
        CAMID=EO;
        VETO_ACTIVATED=false;
        fprintf('Veto setup: [%g,%g], {',DEPTH_LIMITS);
        for i=1:length(LTR)
            fprintf('[');
            fprintf('%d,',LTR{i});
            fprintf(']');
        end
        fprintf('}\n');
        return;
    else
        v=VETO_ACTIVATED;
        return;
    end
end
% Calculate indices for unknown parameters.
base=1;

[ixIO,base]=pindex(nnz(cIO),base);
[ixEO,base]=pindex(nnz(cEO),base);

if (prod(size(cOP))==1)
	cOP=repmat(cOP,size(OP));
elseif (size(cOP,1)<size(OP,1))
	cOP=repmat(cOP,size(OP,1),1);
end
[ixOP,base]=pindex(nnz(cOP),base);

if (size(cEO,1)<size(EO,1))
	cEO(end+1,1)=0;
end

% Copy the current approximations of the unknown values.
IO(cIO)=x(ixIO);
EO(cEO)=x(ixEO);
OP(cOP)=x(ixOP);

v=false;

% Check if any point violated the depth limits.
d=pm_multidepth(IO,EO,OP,vis,1);

bad=d<DEPTH_LIMITS(1) | d>DEPTH_LIMITS(2);
for i=find(any(bad))
    depth=d(:,i);
    j=find(bad(:,i));
    if any(j)
        VETO_ACTIVATED=true;
        v=true;
        return;
        for jj=1:length(j)
            fprintf('Veto: Distance violated for point %d, camera %d: %g (%g,%g)\n',...
                    OPID(j(jj)),CAMID(i),depth(j(jj)),DEPTH_LIMITS);
        end
    end
end

% Check if any cameras are in the wrong order.
if ~isempty(LTR)
    % Calculate center of object.
    c=mean(OP(:,OPIX),2);
    for i=1:length(LTR)
        ltr=LTR{i};
        % Vector for left camera.
        v1=EO(1:3,ltr(1:end-1))-repmat(c,1,length(ltr)-1);
        % Vector for right camera.
        v2=EO(1:3,ltr(2:end))-repmat(c,1,length(ltr)-1);
        % Cross product should point in the same half-plane as UPDIR
        dir=cross(v1,v2)'*UPDIR;
        if any(dir<0)
            VETO_ACTIVATED=true;
            v=true;
            i=find(dir<0,1);
            fprintf('Veto: Order violated for cameras %d and %d\n',...
                    CAMID(ltr(i)),CAMID(ltr(i+1)));
            return;
        end
    end
end
