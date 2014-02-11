function [x,code,n,r,J,T,rr,alphas]=gauss_newton_armijo(...
    resFun,vetoFun,x0,maxIter,convTol,trace,sTest,mu,alphaMin,params)
%GAUSS_NEWTON_ARMIJO Gauss-Newton-Armijo least squares adjustment algorithm.
%
%   [X,CODE,I]=GAUSS_NEWTON_ARMIJO(RES,VETO,X0,N,TOL,TRACE,STEST,MU,AMIN,PARAMS)
%   runs the Gauss-Newton-Armijo least squares adjustment algorithm on the
%   problem with residual function RES and with initial values X0. A maximum
%   of N iterations are allowed and the convergence tolerance is TOL. The
%   final estimate is returned in X. If STEST is true, the iterations are
%   terminated if a (near) singularity warning on the normal matrix is
%   issued. The steplength algorithm uses the constant MU for the Armijo
%   condition and accepts steplengths down to AMIN. In addition, if supplied
%   and non-empty, the VETO function is called to verify that the steplength
%   does not generate an invalid point. The number of iteration I and a
%   success code (0 - OK, -1 - too many iterations, -2 - matrix is singular,
%   -3 - no alpha found) are also returned. If TRACE is true, output sigma0
%   estimates at each iteration.
%
%   [X,CODE,I,R,J]=... also returns the final estimates of the residual
%   vector R and Jacobian matrix J.
%
%   [X,CODE,I,R,J,T,RR,ALPHAS]=... returns the iteration trace as successive
%   columns in T, the successive estimates of sigma0 in RR, and the used
%   steplengths in ALPHAS.
%
%   The function RES is assumed to return the residual function and its
%   Jacobian when called [R,J]=feval(RES,X0,PARAMS{:}), where the cell array
%   PARAMS contains any extra parameters for the residual function.
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%     Nocedal, Wright (2006), "Numerical Optimization", 2nd ed.
%       Springer. ISBN 978-0-387-40065-5.
%     Armijo (1966), "Minimization of Functions Having Lipschitz
%       Continuous First Partial Derivatives", Pacific Journal of
%       Mathematics, 16(1):1-3.
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

% OK until signalled otherwise.
code=0;

% Residual norm trace.
rr=[];

% Steplength vector.
alphas=[];

while true
    % Calculate residual and Jacobian at current point.
    [r,J]=feval(resFun,x,params{:});

    rr(end+1)=sqrt(r'*r);
    if trace
        if isempty(alphas)
            fprintf('Gauss-Newton-Armijo: iteration %d, residual norm=%.2g\n', ...
                    n,rr(end))
        else
            fprintf(['Gauss-Newton-Armijo: iteration %d, residual norm=%.2g, ' ...
                     'last alpha=%s\n'],n,rr(end),...
                    fliplr(deblank(fliplr(deblank(rats(alphas(end)))))));
        end
    end
    
    % Solve normal equations.
    p=(J'*J)\-(J'*r);

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
    % where norm(r) is close to zero at the solution.
    if norm(Jp)<=convTol*norm(r)
        % Converged.
        break
    end

    % Update iteration count.
    n=n+1;
    
    % Perform line search along p.
    [alpha,xNew,rNew]=linesearch(resFun,vetoFun,x,p,alphaMin,r,r'*Jp,mu,params);
	
    % Always update current point and residual. 
    x=xNew;
    r=rNew;

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
        % Ensure that length of residual norm vector matches length of  alpha.
        rr(end+1)=sqrt(r'*r);
        break;
    end

end

% Trim unused trace columns.
if nargout>5
    T=T(:,1:n+1);
end


function [alpha,xNew,rNew]=linesearch(fun,veto,x,p,alphaMin,r0,fp0,mu,params)
%LINESEARCH Perform Armijo linesearch with backtracking.

% Calculate current objective function value.
f0=1/2*(r0'*r0);

% Always try full step length first.
alpha=1;

% Continue while step length is long enough
while alpha>=alphaMin
    % Examine residual at trial point.
    t=x+alpha*p;
    r=feval(fun,t,params{:});
    f=1/2*(r'*r);

    % Is the reduction large enough to satisfy the Armijo condition?
    redOK=f<f0+mu*alpha*fp0;

    if redOK && ~isempty(veto)
        % Call veto function only if reduction is large enough.
        fail=feval(veto,t,params{:});
    else
        fail=false;
    end
    
    if redOK && ~fail
        % If the reduction is large enough and the point did not fail the
        % veto test, we're done.
        xNew=t;
        rNew=r;
        return;
    else
        % Otherwise, try with half the step length.
        alpha=alpha/2;
    end
end

% No acceptable reduction found. Return alpha=0.
alpha=0;
xNew=x;
rNew=f0;
