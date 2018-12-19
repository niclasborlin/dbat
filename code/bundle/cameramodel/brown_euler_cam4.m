function [f,J,JJ]=brown_euler_cam4(x,s)
%BROWN_EULER_CAM4 Residual function for the Brown camera model with Euler angles.
%
%   R=BROWN_EULER_CAM4(X,S) returns the residual vector R of the camera
%   network defined in S evaluated with the approximate values in X.
%
%   Use [R,J]=... or [R,J,JJ]=... to get the analytical Jacobian J
%   and the numerical Jacobian JJ.
%
%   The internal camera model is from Brown (1971) with K1, K2, K3,
%   P1, and P2. The external orientation uses Euler omega-phi-kappa
%   angles.  The values of S.estIO, S.estEO, and S.estOP indicate what
%   parameters in S.IO, S.EO, and S.OP, respectively are unknown.
%
%   References: Brown (1971), "Close-range camera calibration".
%       Photogrammetric Engineering, 37(8): 855-866.
%
%   See also PROB2DBATSTRUCT, BUNDLE, GAUSS_NEWTON_ARMIJO,
%       LEVENBERG_MARQUARDT, LEVENBERG_MARQUARDT_POWELL.

% Update DBAT structure with current estimates in x.
s=deserialize(s,x);

if (nargout>2)
    % Numerical approximation of Jacobian (for debugging only).
    fun=@(x)feval(mfilename,x,s);
    JJ=jacapprox(fun,x);
end

distModel=unique(s.IO.model.distModel);
if ~isscalar(distModel)
    error('Mixed lens distortion models not implemented.');
end

switch distModel
  case 1 % backward/photogrammetric
        
    % Legacy modifications.
    s.bundle.est.IO(:,2:end)=false;
    s.EO.cam(:)=1;
    s.IP.cam(:)=1;
    if (nargout<2)
        % Only residual vector requested.
        
        % Project into pinhole camera.
        xy=multieulerpinhole(s.IO.val,s.IO.model.nK,s.IO.model.nP, ...
                             s.EO.val,s.EO.cam,s.OP.val,s.IP.vis);

        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.IP.val,s.IO.sensor.pxSize,s.IP.cam);
        
        % Compute lens distortion for all measured point.
        ld=multilensdist(m,s.IO.val,s.IO.model.nK,s.IO.model.nP,s.IP.cam);
        
        % Remove lens distortion from measured points.
        ptCorr=m-ld;
        
        % Compute residual for image observations.
        fObs=xy-ptCorr;

        % Residuals for prior observations.
        fPre=prior_obs(x,s,true);
        
        % Pre-allocate residual vector.
        f=nan(s.post.res.ix.n,1);
        % Insert image residuals...
        f(s.post.res.ix.IP)=fObs(:);
        % ...IO residuals...
        f(s.post.res.ix.IO)=fPre.IO;
        % ...EO residuals...
        f(s.post.res.ix.EO)=fPre.EO;
        % ...OP residuals...
        f(s.post.res.ix.OP)=fPre.OP;
    else
        % Project into pinhole camera.
        [xy,dIO1,dEO,dOP]=multieulerpinhole(s.IO.val,s.IO.model.nK, ...
                                            s.IO.model.nP,s.EO.val, ...
                                            s.EO.cam, s.OP.val, ...
                                            s.IP.vis, ...
                                            s.bundle.est.IO, ...
                                            s.bundle.est.EO, s.bundle.est.OP);
        
        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.IP.val,s.IO.sensor.pxSize,s.IP.cam);

        % Compute lens distortion for all measured point.
        [ld,dIO2]=multilensdist(m,s.IO.val,s.IO.model.nK, ...
                                s.IO.model.nP,s.IP.cam,s.bundle.est.IO);

        % Remove lens distortion from measured points.
        ptCorr=m-ld;
        
        % Compute residual for image observations.
        fObs=xy-ptCorr;

        % Compute residual for prior observations.
        [fPre,Jpre]=prior_obs(x,s,true);

        % Pre-allocate residual vector.
        f=nan(s.post.res.ix.n,1);
        % Insert image residuals...
        f(s.post.res.ix.IP)=fObs(:);
        % ...IO residuals...
        f(s.post.res.ix.IO)=fPre.IO;
        % ...EO residuals...
        f(s.post.res.ix.EO)=fPre.EO;
        % ...OP residuals...
        f(s.post.res.ix.OP)=fPre.OP;

        J=sparse(length(f),length(x));
        % Insert image point jacobians...
        J(s.post.res.ix.IP,s.bundle.serial.IO.dest)=dIO1+dIO2;
        J(s.post.res.ix.IP,s.bundle.serial.EO.dest)=dEO;
        J(s.post.res.ix.IP,s.bundle.serial.OP.dest)=dOP;
        % ...IO jacobian...
        J(s.post.res.ix.IO,:)=Jpre.IO;
        % ...EO jacobian...
        J(s.post.res.ix.EO,:)=Jpre.EO;
        % ...OP jacobian...
        J(s.post.res.ix.OP,:)=Jpre.OP;
    end
  case {2,3,4,5}
    % Compute residuals for each camera.
    if (nargout<2)
        % Only residual vector requested.
        
        funs={@res_euler_brown_0,@res_euler_brown_1,...
              @res_euler_brown_2,@res_euler_brown_3};
        % Compute residual for image observations.
        fObs=multi_res(s,funs{distModel-1});

        % Residuals for prior observations.
        fPre=prior_obs(x,s,true);
        
        % Pre-allocate residual vector.
        f=nan(s.post.res.ix.n,1);
        % Insert image residuals...
        f(s.post.res.ix.IP)=fObs(:);
        % ...IO residuals...
        f(s.post.res.ix.IO)=fPre.IO;
        % ...EO residuals...
        f(s.post.res.ix.EO)=fPre.EO;
        % ...OP residuals...
        f(s.post.res.ix.OP)=fPre.OP;
    else
        %fun=@(x)feval(mfilename,x,s);
        % f=feval(fun,x);
        %JJ=jacapprox(fun,x);
        % J=JJ;
        % return;

        funs={@res_euler_brown_0,@res_euler_brown_1,...
              @res_euler_brown_2,@res_euler_brown_3};

        % Compute residual and Jacobian for image
        % observations.
        [fObs,JObs]=multi_res(s,funs{distModel-1});
        
        % Compute residual for prior observations.
        [fPre,Jpre]=prior_obs(x,s,true);

        % Pre-allocate residual vector.
        f=nan(s.post.res.ix.n,1);
        % Insert image residuals...
        f(s.post.res.ix.IP)=fObs(:);
        % ...IO residuals...
        f(s.post.res.ix.IO)=fPre.IO;
        % ...EO residuals...
        f(s.post.res.ix.EO)=fPre.EO;
        % ...OP residuals...
        f(s.post.res.ix.OP)=fPre.OP;

        J=spalloc(length(f),s.bundle.serial.n,nnz(JObs)+nnz(Jpre.IO)+...
                  nnz(Jpre.EO)+nnz(Jpre.OP));
        % Insert image point jacobians...
        J(s.post.res.ix.IP,:)=JObs;
        % ...IO jacobian...
        J(s.post.res.ix.IO,:)=Jpre.IO;
        % ...EO jacobian...
        J(s.post.res.ix.EO,:)=Jpre.EO;
        % ...OP jacobian...
        J(s.post.res.ix.OP,:)=Jpre.OP;
    end
  case -1 % forward/computer vision
    % Legacy modifications.
    s.bundle.est.IO(:,2:end)=false;
    s.EO.cam(:)=1;
    s.IP.cam(:)=1;
    if (nargout<2)
        % Only residual vector requested.
        
        % Project into pinhole camera.
        xy=multieulerpinhole(s.IO.val,s.IO.model.nK,s.IO.model.nP, ...
                             s.EO.val,s.EO.cam,s.OP.val,s.IP.vis);

        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.IP.val,s.IO.sensor.pxSize,s.IP.cam);
        
        % Compute lens distortion for projected points.
        ld=multilensdist(xy,s.IO.val,s.IO.model.nK,s.IO.model.nP,s.IP.cam);

        % Add lens distortion to projected points.
        ptDist=xy+ld;
        
        % Compute residual for prior observations.
        fPre=pm_preobs(x,s);
    
        f=[ptDist(:)-m(:);fPre];
    else
        % Project into pinhole camera.
        [xy,dIO1,dEO,dOP]=multieulerpinhole(s.IO.val,s.IO.model.nK, ...
                                            s.IO.model.nP,s.EO.val, ...
                                            s.EO.cam, s.OP.val, ...
                                            s.IP.vis, ...
                                            s.bundle.est.IO, ...
                                            s.bundle.est.EO,s.bundle.est.OP);
        
        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.IP.val,s.IO.sensor.pxSize,s.IP.cam);

        % Create arrays of columns indices for IO derivatives.
        [ixpp,ixf,ixK1,ixP1]= ...
            createiocolumnindices(s.bundle.est.IO(:,1), ...
                                  s.IO.model.nK,s.IO.model.nP);
        
        % No need to compute the partials w.r.t. pp or f.
        est=s.bundle.est.IO(:,1);
        if nnz(ixpp)
            est(ixpp)=0;
        end
        if nnz(ixf)
            est(ixf)=0;
        end
        % Potentially new column indices for lens distortion parameters.
        [~,~,ixK2,ixP2]=createiocolumnindices(est,s.IO.model.nK,s.IO.model.nP);

        % Compute lens distortion for projected points.
        [ld,dIO2,dxy]=multilensdist(xy,s.IO.val(:,1),s.IO.model.nK, ...
                                    s.IO.model.nP,s.IP.cam,est);
        
        % Add lens distortion to projected points.
        ptDist=xy+ld;
        
        % Compute residual for prior observations.
        [fPre,Jpre]=prior_obs(x,s,true);

        fObs=ptDist(:)-m(:);

        % Pre-allocate residual vector.
        f=nan(s.post.res.ix.n,1);
        % Insert image residuals...
        f(s.post.res.ix.IP)=fObs(:);
        % ...IO residuals...
        f(s.post.res.ix.IO)=fPre.IO;
        % ...EO residuals...
        f(s.post.res.ix.EO)=fPre.EO;
        % ...OP residuals...
        f(s.post.res.ix.OP)=fPre.OP;

        % Combine all lens distortion parameters
        ixLD1=[ixK1;ixP1];
        ixLD2=[ixK2;ixP2];
        if nnz(ixLD1)>0
            % Insert partials w.r.t. K and P into dIO1.
            dIO1(:,ixLD1(ixLD1~=0))=dIO2(:,ixLD2(ixLD2~=0));
        end
        
        if nnz(ixf)>0
            % Update Jacobian w.r.t. focal length.
            dIO1(:,ixf)=dIO1(:,ixf)+dxy*dIO1(:,ixf);
        end
        
        % Update Jacobian w.r.t. EO and OP.
        dEO2=dEO+dxy*dEO;
        dOP2=dOP+dxy*dOP;
        
        J=sparse(length(f),length(x));
        % Insert image point jacobians...
        J(s.post.res.ix.IP,s.bundle.serial.IO.dest)=dIO1;
        J(s.post.res.ix.IP,s.bundle.serial.EO.dest)=dEO2;
        J(s.post.res.ix.IP,s.bundle.serial.OP.dest)=dOP2;
        % ...IO jacobian...
        J(s.post.res.ix.IO,:)=Jpre.IO;
        % ...EO jacobian...
        J(s.post.res.ix.EO,:)=Jpre.EO;
        % ...OP jacobian...
        J(s.post.res.ix.OP,:)=Jpre.OP;
    end
  otherwise
    error('Bad distortion model %d',distModel);
end

