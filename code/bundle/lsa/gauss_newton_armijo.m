function [x,code,n,f,J,T,rr,alphas]=gauss_newton_armijo(...
    resFun,vetoFun,x0,mu,alphaMin, maxIter,convTol,trace, sTest,params)
%GAUSS_NEWTON_ARMIJO Gauss-Newton-Armijo least squares adjustment algorithm.
%
%   [X,CODE,I]=GAUSS_NEWTON_ARMIJO(RES,VETO,X0,MU,AMIN,N,TOL,TRACE,STEST,PARAMS)
%   runs the Gauss-Newton-Armijo least squares adjustment algorithm on the
%   problem with residual function RES and with initial values X0. A maximum
%   of N iterations are allowed and the convergence tolerance is TOL. The
%   final estimate is returned in X. If STEST is true, the iterations are
%   terminated if a (near) singularity warning on the normal matrix is
%   issued. The steplength algorithm uses the constant MU for the Armijo
%   condition and accepts steplengths down to AMIN. In addition, the VETO
%   function is called to verify that the steplength does not generate an
%   invalid point. The number of iteration I and a success code (0 - OK, -1
%   - too many iterations, -2 - matrix is singular, -3 - no alpha found) are
%   also returned. If TRACE is true, output sigma0 estimates at each
%   iteration.
%
%   [X,CODE,I,F,J]=... also returns the final estimates of the residual
%   vector F and jacobian matrix J.
%
%   [X,CODE,I,F,J,T,RR,ALPHAS]=... returns the iteration trace as successive
%   columns in T, the successive estimates of sigma0 in RR, and the used
%   steplengths in ALPHAS.
%
%   The function RES is assumed to return the residual function and its
%   jacobian when called [F,J]=feval(RES,X0,PARAMS{:}), where the cell array
%   PARAMS contains any extra parameters for the residual function.
%
%See also: BUNDLE, GAUSS_MARKOV, LEVENBERG_MARQUARDT,
%   LEVENBERG_MARQUARDT_POWELL.

% $Id$

% Initialize current estimate and iteration trace.
x=x0;

if nargout>5
    % Pre-allocate fixed block if trace is asked for.
    blockSize=50;
    T=nan(length(x),min(blockSize,maxIter+1));
    % Enter x0 as first column.
    T(:,1)=x0;
end

if sTest
    % Clear last warning.
    lastwarn('');
end

% Iteration counter.
n=0;

% OK until proven otherwise.
code=0;

% Residual norm trace.
rr=[];

% Steplength vector.
alphas=[];

while true
    % Calculate residual and jacobian at current point.
    [f,J]=feval(resFun,x,params{:});

    rr(end+1)=sqrt(f'*f);
    if trace
        if isempty(alphas)
            fprintf('Gauss-Newton-Armijo: iteration %d, residual norm=%.1g\n', ...
                    n,rr(end))
        else
            fprintf(['Gauss-Newton-Armijo: iteration %d, residual norm=%.1g, ' ...
                     'last alpha=%s\n'],n,rr(end),...
                    fliplr(deblank(fliplr(deblank(rats(alphas(end)))))));
        end
    end
    
    % Solve normal equations.
    p=(J'*J)\-(J'*f);

    if sTest
        % Check if we got the 'Matrix is singular to working precision' warning.
        [warnmsg,msgid]=lastwarn; %#ok<ASGLU>
        if strcmp(msgid,'MATLAB:singularMatrix') || ...
                strcmp(msgid,'MATLAB:nearlySingularMatrix')
            code=-2;
            break
        end
    end

    % Pre-calculate J*p.
    Jp=J*p;
    
    % Terminate if angle between projected residual is smaller than
    % threshold. Warning! This test may be very strict on synthetic data
    % where norm(f) is close to zero at the solution.
    if norm(Jp)<=convTol*norm(f)
        % Converged.
        break
    end

    % Update iteration count.
    n=n+1;
    
    % Perform line search along p.
    [alpha,xNew,fNew]=linesearch(resFun,vetoFun,x,p,alphaMin,f,f'*Jp,mu,params);
	
    % Always update current point and residual. 
    x=xNew;
    f=fNew;

    % Store iteration trace and algorithm performance parameters.
    alphas(end+1)=alpha; %#ok<AGROW>
    if nargout>5
        % Store iteration trace.
        if n+1>size(T,2)
            % Expand by blocksize if needed.
            T=[T,nan(length(x),blockSize)]; %#ok<AGROW>
        end
        T(:,n+1)=x;
    end

    if alpha==0
        % Not sufficient descent, stop.
        code=-2;
        % Ensure that length of residual norm vector matches length of alpha.
        rr(end+1)=rr(end);
        break;
    end
    
    if n>maxIter
        % Too many iterations, stop.
        code=-1; %
        % Ensure that length of residual norm vector matches length of alpha.
        rr(end+1)=nan;
        break;
    end

end

% Trim unused trace columns.
if nargout>5
    T=T(:,1:n+1);
end


function [alpha,xNew,fNew]=linesearch(fun,veto,x,p,alphaMin,f0,fp0,mu,params)
%LINESEARCH Perform Armijo linesearch with backtracking.

% Calculate current objective function value.
obj0=0.5*(f0'*f0);

% Always start with full step length.
alpha=1;

% Continue while step length is long enough
while alpha>=alphaMin
    % Examine residual at trial point.
    trial=x+alpha*p;
    res=feval(fun,trial,params{:});
    obj=0.5*(res'*res);

    % Is the reduction large enough to satisfy the Armijo condition?
    redOK=obj<obj0+mu*alpha*fp0;

    if redOK && ~isempty(veto)
        % Call veto function only if reduction is large enough.
        fail=feval(veto,trial,params{:});
    else
        fail=false;
    end
    
    if redOK && ~fail
        % If the reduction is large enough and the point did not fail the veto
        % function, we're done.
        xNew=trial;
        fNew=obj;        
        return;
    else
        % Otherwise, try with half the step length.
        alpha=alpha/2;
    end
end

% No acceptable reduction found. Return alpha=0.
alpha=0;
xNew=x;
fNew=obj0;
