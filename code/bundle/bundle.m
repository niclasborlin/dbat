function [s,ok,iters,s0,E]=bundle(s,varargin)
%BUNDLE Run bundle adjustment iterations on a camera network.
%
%   [S,OK,N]=BUNDLE(S), where S is a struct returned by PROB2DBATSTRUCT,
%   runs the damped bundle adjustment on the camera network in the structure
%   S. The parameter values in S are used as initial values. The cIO, cEO,
%   cOP fields of S are used to indicate which parameters are free. OK is
%   returned as true if the bundle converged within the allowed number of
%   iterations. N gives the number of iterations. On return, the parameter
%   values in S are updated if the bundle converged.
%
%   ...=BUNDLE(S,...,K), where K is an integer, sets the maximum number of
%   iterations to K (default: 20).
%
%   ...=BUNDLE(S,...,DAMP), where DAMP is a string, specifies which damping
%   to use: 'none' or 'GM' (classic bundle with no damping), 'GNA'
%   (Gauss-Newton with Armijo linesearch, default), 'LM' (original
%   Levenberg-Marquardt) , 'LMP' (Levenberg-Marquardt with Powell dogleg).
%
%   ...=BUNDLE(S,...,'trace') specifies that the bundle should
%   print an trace during the iterations.
%
%   ...=BUNDLE(S,...,TOL), where TOL is a scalar < 1, sets the
%   convergence tolerance. Default: 1e-06.
%
%   ...=BUNDLE(S,...,'absterm') specifies that the bundle should
%   use an absolute termination criteria instead of a relative (the
%   default). This might be useful if you are testing the bundle on
%   synthetic data with zero residual.
%
%   ...=BUNDLE(S,...,'singulartest') specifies that the bundle should
%   stop immediately if a 'Matrix is singular' or 'Matrix is almost
%   singular' warning is issued on the normal matrix. This is the
%   default. Use ...=BUNDLE(S,...,'nosingulartest') to inhibit
%   singular tests.
%
%   ...=BUNDLE(S,...,CHI), where CHI is a logical scalar, specifies if
%   chirality veto damping should be used (default: false). Chirality
%   veto damping is ignored for the undamped bundle.
%
%   ...=BUNDLE(S,...,'pmdof') uses Photomodeler degrees-of-freedom
%   computation.
%
%   ...=BUNDLE(S,...,'dofverb') outputs how the degrees of freedom
%   are calculated.
%
%   The used distortion model is specified by the vector
%   s.IOdistModel. The available values are:
%     1 - Legacy Photogrammetry, no affine (slightly faster than
%         2+). Default for DBAT versions before v0.7.1.
%     2 - Flexible Photogrammetry, no affine (replica of 1)
%     3 - Photogrammetry, affine before lens distortion. Default
%         since DBAT versions v0.7.2. 
%     4 - Photogrammetry, affine after lens distortion
%     5 - Photogrammetry, anisotropic scale before lens dist, skew after.
%    -1 - Computer Vision, no affine (fast).
%
%   [S,OK,N,S0]=... returns the sigma0 (in pixel units) for the last iteration.
%
%   [S,OK,N,S0,E]=... returns a struct E with information about the iterations:
%       E.trace   - NOBS-by-(N+1) array with successive parameter estimates.
%       E.res     - (N+1)-vector with the residual norm at every iteration.
%       E.damping - struct with damping-specific information, including
%           E.damping.name - the name of the damping scheme used.
%
%   Use BUNDLE_COV to compute covariances, etc., of the result or
%   BUNDLE_RESULT_FILE to generate a result file.
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%
%See also: GAUSS_MARKOV, GAUSS_NEWTON_ARMIJO, LEVENBERG_MARQUARDT,
%   LEVENBERG_MARQUARDT_POWELL, BROWN_EULER_CAM2, BUNDLE_COV,
%   PROB2DBATSTRUCT, BUNDLE_RESULT_FILE

maxIter=20;
damping='gna';
veto=false;
singularTest=true;
doTrace=false;
dofVerb=false;
pmDof=false;
absTerm=false;
convTol=1e-6;

while ~isempty(varargin)
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        % N
        v=varargin{1};
        if v==round(v)
            maxIter=v;
        else
            convTol=v;
        end
        varargin(1)=[];
    elseif ischar(varargin{1})
        % DAMP, TRACE, or CXX
        switch lower(varargin{1})
          case {'none','gm','gna','lm','lmp'}
            % OK
            damping=varargin{1};
            varargin(1)=[];
          case 'trace'
            doTrace=true;
            varargin(1)=[];
          case 'singulartest'
            singularTest=true;
            varargin(1)=[];
          case 'nosingulartest'
            singularTest=false;
            varargin(1)=[];
          case 'pmdof'
            pmDof=true;
            varargin(1)=[];
          case 'dofverb'
            dofVerb=true;
            varargin(1)=[];
          case 'absterm'
            absTerm=true;
            varargin(1)=[];
          otherwise
            error('DBAT:bundle:badInput','Unknown damping');
        end
    elseif islogical(varargin{1})
        veto=varargin{1};
        varargin(1)=[];
    else
        error('DBAT:bundle:badInput','Unknown parameter');
    end
end

% Verify a consistent estimation/use observations arrays. If we
% shouldn't estimate a parameter, we should not use the prior
% observation of it during the bundle.
if any(s.prior.IO.use(~s.bundle.est.IO))
    warning('DBAT:bundle:badInput',...
            'Some IO parameters are set to both ''fixed'' and ''observed''');
    disp('Setting IO parameters to fixed');
    s.prior.IO.use(~s.bundle.est.IO)=false;
end
if any(s.prior.EO.use(~s.bundle.est.EO))
    warning('DBAT:bundle:badInput',...
            'Some EO parameters are set to both ''fixed'' and ''observed''');
    disp('Setting EO parameters to fixed');
    s.prior.EO.use(~s.bundle.est.EO)=false;
end
if any(s.prior.OP.use(~s.bundle.est.OP))
    warning('DBAT:bundle:badInput',...
            'Some OP parameters are set to both ''fixed'' and ''observed''');
    disp('Setting OP parameters to fixed');
    s.prior.OP.use(~s.bundle.est.OP)=false;
end

if isempty(s.bundle.serial) || isempty(s.bundle.deserial)
    % Compute (de-)serialization indices.
    s=buildserialindices(s);
end

% Create x0 vector.
[x0,paramTypes]=serialize(s);

% Residual function.
resFun=@(x)brown_euler_cam4(x,s);
%resFun=@(x)both_brown_res(x,s);

if veto
    vetoFun=@chirality;
else
    vetoFun='';
end

% Compute weight matrix.
W=buildweightmatrix(s);

% Choose between a relative and absolute termination criteria. The
% relative termination criteria (default) corresponds to the angle
% between the residual vector and the tangent plane of the non-linear
% residual surface. For real-world, noisy data, it is generally the
% best, as it does not depend on the scaling of the parameters or
% observations. However, for synthetic data with zero residual, the
% termination critera may never activate and the bundle will fail. In
% such cases, an absolute termination criteria on the residual vector
% alone can be used.
if absTerm
    % Absolute termination criteria.
    termFun=@(Jp,r)norm(r)<=convTol;
else
    % Relative termination criteria (angle).
    termFun=@(Jp,r)norm(Jp)<=convTol*norm(r);
end

% For all optimization methods below, the final estimate is returned
% in x.  The final weighted and unweighted residual vectors and
% Jacobians are returned in the struct final. Successive estimates of x and
% norm(r) are returned as columns of X and elements of res,
% respectively. Furthermore, a status code (0 - ok, -1 - too many
% iterations) and the number of required iterations are returned.

% Set up return struct with bundle setup.

% TODO: Verify consistent return sizes of trace array, residual
% estimates, and damping et al. vectors given
% 1) x* = x0
% 2) x* found before maxIter
% 3) x* found at maxIter
% 4) x* not found at maxIter

distModel=unique(s.IO.model.distModel);

% Display lens distortion models.
if isscalar(distModel) && distModel>0
    fprintf(['Using Backward Brown (Photogrammetry) lens distortion ' ...
             'model %d for all cameras\n'],distModel);
elseif isscalar(distModel) && distModel<0
    fprintf(['Using Forward Brown (Computer Vision) lens distortion ' ...
             'model %d for all cameras\n'],-distModel);
else
    disp('Using mix of Forward/Backward Brown lens distortion models');
end

% Display any self-calbration parameters.
selfCal=any(s.bundle.est.IO,1);
strs=s.IO.type(:,1);
if all(selfCal)
    allParamCal=all(s.bundle.est.IO,2);
    anyParamCal=any(s.bundle.est.IO,2);
    if all(allParamCal==anyParamCal)
        % Any parameters that is estimated in one camera is
        % estimated in all cameras.
        calParams=sprintf('%s ',strs{allParamCal(1:length(strs))});
        disp(['Self-calibration: yes (',strtrim(calParams),')']);
    else
        disp('Self-calibration: yes (mixed parameters)');
    end
elseif ~any(selfCal)
    disp('Self-calibration: no');
else
    disp('Self-calibration: mixed');
    warning('Mixed self-calibration is poorly tested.')
end

% Warn if non-zero aspect/skew is specified for a model that does
% not support it.
modelsWithoutB=ismember(s.IO.model.distModel,1:2);
usedBadModels=unique(s.IO.model.distModel(modelsWithoutB));
if any(s.IO.val(4:5,modelsWithoutB)~=0)
    warning(['Non-zero aspect and/or skew specified. This is not ' ...
             'supported by lens distortion model %d! Results may ' ...
             'be inaccurate.'],usedBadModels);
end
% Warn if asked to estimate aspect and/or skew with a model that
% does not suppos.Irt it.
if any(s.bundle.est.IO(4:5,modelsWithoutB))
    warning(['Trying to estimate aspect and/or skew. This is not ' ...
             'supported by lens distortion model %d! Results will ' ...
             'be inaccurate.'],usedBadModels);
end
    
% Version string.
[v,d]=dbatversion;
E=struct('maxIter',maxIter,'convTol',convTol,'absTerm',absTerm, ...
         'singularTest',singularTest, 'chirality',veto,'dateStamp', ...
         datestr(now), 'version',sprintf('%s (%s)',v,d));

switch lower(damping)
  case {'none','gm'}
    % Gauss-Markov with no damping.
    
    % Call Gauss-Markov optimization routine.
    stopWatch=cputime;
    [x,code,iters,final,X,res]=gauss_markov(resFun,x0,W,maxIter, ...
                                          termFun,doTrace, singularTest);
    time=cputime-stopWatch;
    E.damping=struct('name','gm');
  case 'gna'
    % Gauss-Newton with Armijo linesearch.

    % Armijo parameter.
    mu=0.1;
    % Shortest allowed step length.
    alphaMin=1e-9;
    
    % Call Gauss-Newton-Armijo optimization routine. The vector alpha is
    % returned with the step lengths used at each iteration.
    stopWatch=cputime;
    [x,code,iters,final,X,res,alpha]=gauss_newton_armijo(resFun, ...
                                                      vetoFun,x0,W, ...
                                                      maxIter, ...
                                                      termFun,doTrace, ...
                                                      singularTest, ...
                                                      mu,alphaMin);  
    time=cputime-stopWatch;
    E.damping=struct('name','gna','alpha',alpha,'mu',mu,'alphaMin',alphaMin);
  case 'lm'
    % Original Levenberg-Marquardt "lambda"-version.

    % Initial lambda value. A negative value tells levenberg_marquardt to
    % scale lambda0 by trace(J0'*J0)/size(J0,2).
    lambda0=-1e-10;

    % Use lambda0 as lower damping cutoff value.
    lambdaMin=lambda0;

    % Call Levenberg-Marquardt optimization routine. The vector lambdas is
    % returned with the lambda values used at each iteration.
    stopWatch=cputime;
    [x,code,iters,final,X,res,lambda]=levenberg_marquardt(resFun,vetoFun,x0, ...
                                                      W,maxIter, ...
                                                      termFun,doTrace, ...
                                                      lambda0,lambdaMin);
    time=cputime-stopWatch;
    E.damping=struct('name','lm','lambda',lambda,'lambda0',lambda(1), ...
                     'lambdaMin',lambda(1));
  case 'lmp'
    % Levenberg-Marquardt-Powell trust-region, "delta"-version with
    % Powell dogleg.
    
    % Limits on the gain ratio.
    rhoBad=0.25; % Below this is bad.
    rhoGood=0.75; % Above this is good.
    
    % Initial delta value.
    delta0=norm(x0);
    
    % Call Levenberg-Marquardt-Powell optimization routine. The vector deltas
    % and rhos are returned with the used delta and computed rho values for
    % each iteration.
    stopWatch=cputime;
    [x,code,iters,final,X,res,delta,rho,step]=levenberg_marquardt_powell(...
        resFun,vetoFun,x0,W,maxIter,termFun,doTrace, delta0,rhoBad,rhoGood);
    time=cputime-stopWatch;
    E.damping=struct('name','lmp','delta',delta,'rho',rho,'delta0',delta0,...
                     'rhoBad',rhoBad,'rhoGood',rhoGood,'step',step);
  otherwise
    error('DBAT:bundle:internal','Unknown damping');
end

% Store iteration results.
E.res=res;
E.trace=X;
E.time=time;
E.code=code;
E.usedIters=iters;

% Store final weighted residual and Jacobian for later covariance calculations.
E.final=final;
% Signal that posterior covariance matrix has not been factorized.
E.final.factorized=[];

% Handle returned values.
ok=code==0;

% Update s if optimization converged.
if ok
    s=deserialize(s,x);
end
   
% Update formats.
aspect=ones(2,size(s.IO.val,2));
aspect(1,:)=1+s.IO.val(4,:);
s.post.sensor.imSize=s.IO.sensor.imSize;
s.post.sensor.pxSize=s.IO.sensor.pxSize.*aspect;
imFormat=s.post.sensor.imSize.*s.IO.sensor.pxSize.*aspect;
s.post.sensor.ssSize=imFormat;

E.paramTypes=paramTypes;

% Analyse potential problem with the Jacobian (design matrix).

E.weakness=struct('structural',[],'numerical',[]);

if code==-2
    % Numerically rank deficient, try to figure out why.
    
    E.weakness.numerical.rank=size(E.final.weighted.J,2);
    try
        % Use spnrank to find rank.
        fprintf('Trying to estimate numerical rank of Jacobian...');
        E.weakness.numerical.rank=spnrank(E.final.scaled.J);
        fprintf('done.\n');
    catch
        % If spnrank failed
        E.weakness.numerical.rank=nan;
    end
    % Deficiency.
    E.weakness.numerical.deficiency=...
        size(E.final.scaled.J,2)-E.weakness.numerical.rank;
    if E.weakness.numerical.deficiency>0
        % Find vectors in offending null-space.
        try
            fprintf('Trying to estimate null-space...');
            JTJ=E.final.scaled.J'*E.final.scaled.J;
            % Add shift of sqrt(eps) to diagonal to avoid failure for exact
            % rank-deficiency.
            shift=sqrt(eps);
            JTJ=JTJ+shift*speye(size(JTJ));
            [V,D]=eigs(JTJ,E.weakness.numerical.deficiency,'SM',...
                       struct('issym',true,'isreal',true));
            % Remove shift.
            D=D-shift*speye(size(D));
            fprintf('done.\n');
            % Sort increasingly by eigenvalue.
            [~,i]=sort(abs(diag(D)),'ascend');
            D=D(i,i);
            V=V(:,i);
            E.weakness.numerical.V=V;
            E.weakness.numerical.d=diag(D);
            E.weakness.numerical.trace=trace(JTJ);
            E.weakness.numerical.suspectedParams=cell(1,size(V,2));
            for j=1:size(V,2)
                % Sort values descendingly.
                [~,k]=sort(abs(V(:,j)),'descend');
                v=V(k,j);
                % Keep element halfway between average and max.
                avg=sqrt(1/size(V,1));
                keep=abs(v)>mean([avg,abs(v(1))]);
                k=k(keep);
                v=v(keep);
                sp=struct('values',v,'indices',k,'params',{E.paramTypes(k)});
                E.weakness.numerical.suspectedParams{j}=sp;
            end
        catch
            % If eigs fails.
            E.weakness.numerical.suspectedParams={};
        end
    end
else
    E.weakness.numerical.rank=size(E.final.weighted.J,2);
    E.weakness.numerical.deficiency=0;
end

if code==-4
    % Structural rank deficiency. Record potential cause.
    E.weakness.structural.dmperm=dmperm(E.final.weighted.J);
    E.weakness.structural.rank=nnz(E.weakness.structural.dmperm);
    E.weakness.structural.deficiency=...
        size(E.final.weighted.J,2)-E.weakness.structural.rank;
    E.weakness.structural.suspectedParams=...
        E.paramTypes(E.weakness.structural.dmperm==0);
    
    % Mark numerical rank as unchecked.
    E.weakness.numerical.rank=nan;
    E.weakness.numerical.deficiency=nan;
end

% Always update the residuals.
s.post.res.IP=nan(size(s.IP.val));
s.post.res.IO=nan(size(s.IO.val));
s.post.res.EO=nan(size(s.EO.val));
s.post.res.OP=nan(size(s.OP.val));

s.post.res.IP(:)=final.unweighted.r(s.post.res.ix.IP);
% Mark pt residuals are in mm, scale to pixels.
ptCols=s.IP.ix(s.IP.vis);
s.post.res.IP=s.post.res.IP./s.IO.sensor.pxSize(:,s.IP.cam(ptCols));

% This might be problematic if a block parameter is used as an observation...
s.post.res.IO(s.prior.IO.use)=final.unweighted.r(s.post.res.ix.IO);
s.post.res.EO(s.prior.EO.use)=final.unweighted.r(s.post.res.ix.EO);
s.post.res.OP(s.prior.OP.use)=final.unweighted.r(s.post.res.ix.OP);

% Sigma0 is sqrt(r'*r/(m-n)), where m is the number of
% observations, and n is the number of unknowns.

% Extra observations without a specific reisdual (fixed control points
% that have been measured, fixed camera stations that have been
% used).
if pmDof
    p=nnz(s.bundle.est.OP(:,any(s.IP.vis,2))==0)+nnz(s.bundle.est.EO(1:6,any(s.IP.vis,1))==0);
else
    p=0;
end

r=E.final.weighted.r;
lenR=length(r);
lenX=length(x);
dof=(lenR+p-lenX);
if dofVerb
    fprintf('%s: dof=%d+%d-%d=%d.\n',mfilename,lenR,p,lenX,dof);
end
s0=sqrt((r'*r)/dof);
sigmas=s0*s.IP.sigmas;
s.post.sigmas=sigmas;

E.numObs=lenR;
E.numParams=lenX;
E.redundancy=dof;
E.s0=s0;
E.sigmas=sigmas;


function [f,J,JJ]=both_brown_res(x,s)

if nargout<2
    f=brown_euler_cam(x,s);    
    f2=brown_euler_cam2(x,s);
    v=max(abs(f-f2));
    if v>1e-14
        warning('Residual difference=%g',v);
    end
elseif nargout==2
    [f,J]=brown_euler_cam(x,s);    
    [f2,J2]=brown_euler_cam2(x,s);
    v=max(abs(f-f2));
    if v>1e-14
        warning('Residual difference=%g',v);
    end
    v=full(max(max(abs(J-J2))));
    if v>1e-8
        warning('Jacobian difference=%g',v);
    end
else
    [f,J,JJ]=brown_euler_cam(x,s);
    [f2,J2,JJ2]=brown_euler_cam2(x,s);
    v=max(abs(f-f2));
    if v>1e-14
        warning('Residual difference=%g',v);
    end
    v=full(max(max(abs(J-J2))));
    if v>1e-8
        warning('Jacobian difference=%g',v);
    end
    v=full(max(max(abs(JJ-JJ2))));
    if v>1e-14
        warning('Numerical jacobian difference=%g',v);
    end
end

    
