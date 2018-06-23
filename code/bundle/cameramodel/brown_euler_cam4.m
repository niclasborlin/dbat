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

distModel=unique(s.IOdistModel);
if ~isscalar(distModel)
    error('Mixed lens distortion models not implemented.');
end

if distModel>0
    if distModel==1 % backward/photogrammetric
        
        % Legacy modification.
        s.estIO(:,2:end)=false;
        if (nargout<2)
            % Only residual vector requested.
    
            % Project into pinhole camera.
            xy=multieulerpinhole(s.IO,s.nK,s.nP,s.EO,s.imCams,s.OP,s.vis);

            % Convert measured points from pixels to mm and flip y coordinate.
            m=diag([1,-1])*multiscalepts(s.markPts,s.IO,s.nK,s.nP,s.ptCams);
        
            % Compute lens distortion for all measured point.
            ld=multilensdist(m,s.IO,s.nK,s.nP,s.ptCams);
        
            % Remove lens distortion from measured points.
            ptCorr=m-ld;
            
            % Compute residual for image observations.
            fObs=xy-ptCorr;

            % Compute residual for prior observations.
            fPre=pm_preobs(x,s);

            % Combine residuals.
            f=[fObs(:);fPre(:)];
        else
            % Project into pinhole camera.
            [xy,dIO1,dEO,dOP]=multieulerpinhole(s.IO,s.nK,s.nP,s.EO,s.imCams,...
                                                s.OP,s.vis,s.estIO,s.estEO,...
                                                s.estOP);
            
            % Convert measured points from pixels to mm and flip y coordinate.
            m=diag([1,-1])*multiscalepts(s.markPts,s.IO,s.nK,s.nP,s.ptCams);

            % Compute lens distortion for all measured point.
            [ld,dIO2]=multilensdist(m,s.IO,s.nK,s.nP,s.ptCams,s.estIO);

            % Remove lens distortion from measured points.
            ptCorr=m-ld;
        
            % Compute residual for image observations.
            fObs=xy-ptCorr;

            % Compute residual for prior observations.
            [fPre,Jpre]=pm_preobs(x,s);

            f=[fObs(:);fPre(:)];
            
            J=[dIO1+dIO2,dEO,dOP;Jpre];
        end
    else
        if distModel>=2 && distModel<=5
            % Compute residuals for each camera.
            if (nargout<2)
                % Only residual vector requested.
    
                funs={@res_euler_brown_0,@res_euler_brown_1,...
                      @res_euler_brown_2,@res_euler_brown_3};
                % Compute residual for image observations.
                fObs=multi_res(s.IO,s.EO,s.OP,s,funs{distModel-1});

                % Compute residual for prior observations.
                fPre=pm_preobs(x,s);

                % Combine.
                f=[fObs(:);fPre(:)];
            else
                fun=@(x)feval(mfilename,x,s);
                %JJ=jacapprox(fun,x);

                funs={@res_euler_brown_0,@res_euler_brown_1,...
                      @res_euler_brown_2,@res_euler_brown_3};
                % Compute residual and Jacobian for image
                % observations.
                [fObs,JObs]=multi_res(s.IO,s.EO,s.OP,s,funs{distModel-1});
                
                % Compute residual for prior observations.
                [fPre,Jpre]=pm_preobs(x,s);

                % Combine residual and Jacobian
                f=[fObs(:);fPre(:)];
            
                J=[JObs;Jpre];
            end
        end
    end
elseif all(s.IOdistModel==-1) % forward/computer vision
    if (nargout<2)
        % Only residual vector requested.
    
        % Project into pinhole camera.
        xy=multieulerpinhole(s.IO,s.nK,s.nP,s.EO,s.cams,s.OP,s.vis);

        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.markPts,s.IO,s.nK,s.nP,s.ptCams);
        
        % Compute lens distortion for projected points.
        ld=multilensdist(xy,s.IO,s.nK,s.nP,s.ptCams);

        % Add lens distortion to projected points.
        ptDist=xy+ld;
        
        % Compute residual for prior observations.
        fPre=pm_preobs(x,s);
    
        f=[ptDist(:)-m(:);fPre];
    else
        % Project into pinhole camera.
        [xy,dIO1,dEO,dOP]=multieulerpinhole(s.IO,s.nK,s.nP,s.EO,s.cams,s.OP, ...
                                            s.vis,s.estIO,s.estEO,s.estOP);
        
        % Convert measured points from pixels to mm and flip y coordinate.
        m=diag([1,-1])*multiscalepts(s.markPts,s.IO,s.nK,s.nP,s.ptCams);

        % Create arrays of columns indices for IO derivatives.
        [ixpp,ixf,ixK1,ixP1]=createiocolumnindices(s.estIO,s.nK,s.nP);
        % No need to compute the partials w.r.t. pp or f.
        est=s.estIO;
        if nnz(ixpp)
            est(ixpp)=0;
        end
        if nnz(ixf)
            est(ixf)=0;
        end
        % Potentially new column indices for lens distortion parameters.
        [~,~,ixK2,ixP2]=createiocolumnindices(est,s.nK,s.nP);

        % Compute lens distortion for projected points.
        [ld,dIO2,dxy]=multilensdist(xy,s.IO,s.nK,s.nP,s.ptCams,est);
        
        % Add lens distortion to projected points.
        ptDist=xy+ld;
        
        % Compute residual for prior observations.
        [fPre,Jpre]=pm_preobs(x,s);

        f=[ptDist(:)-m(:);fPre];

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
        
        J=[dIO1,dEO2,dOP2;Jpre];
    end
else
    error('Mixed lens distortion models not implemented.');
end

