function [x,code,n,X,alphas,F]=gaussn_paper(prob,veto,x0,tol,maxIter,alphaMin,mu,W,params)
%GAUSSN Damped Gauss-Newton method, Goldstein-Armijo linesearch.
%
%[x,code,n,X,alphas]=gaussn(prob,veto,x0,tol,maxIter,alphaMin,mu,W,params)
%prob     - name of problem.
%veto     - name of veto function.
%x0       - starting approximation.
%tol      - convergence tolerance.
%maxiter  - maximum number of iterations.
%alphaMin - shortest accepted step length.
%mu       - constant in Armijo's condition.
%W        - weight matrix. Use W=1 for unweighted problems.
%params   - cell array with additional arguments to residual and jacobian funs.
%x        - optimal solution.
%code     - error code
%              0 - OK
%             -1 - Too many iterations
%             -2 - No acceptable reduction found in line search.
%             -3 - Matrix is singular or close to singular.
%n        - number of consumed iterations.
%X        - iteration trace. X(:,i+1) is the ith iteration.
%alphas   - vector of accepted step lengths. alphas(i) is the ith step length.


x=x0;
X=x;
alphas=[];

n=0;

% OK until proven otherwise.
code=0;

% Clear warning.
lastwarn('');

alpha=1;

while true
    % Calculate residual and jacobian at current point.
    [F,J]=feval(prob,x,params{:});
    % Gradient.
    g=J'*F;
    % Calculate Gauss-Newton search direction.
	JTJ=J'*J;
    p=JTJ\-(J'*F);

    % Check if we got the 'Matrix is singular to working precision.' warning.
    [warnmsg,msgid]=lastwarn;
    if strcmp(msgid,'MATLAB:singularMatrix') || ...
            strcmp(msgid,'MATLAB:nearlySingularMatrix')
        code=-3;
        break
    end
    
    % Terminate if angle between projected residual is smaller than
    % threshold. 
    % A cleaner test would be the strict relative test
    % norm(J*p)>tol*norm(F), but that wouldn't work on synthetic data
    % with norm(F) very small at solution.
	%[norm(p),max(abs(p(6:end)))]
	%if (max(abs(p([1,2,6:end])))<tol)
    Jp=J*p;
	if norm(Jp)<tol*norm(F) && alpha==1
        % Converged.
		%X=[X,x+p];
        break;
    end
    n=n+1;
    if n>maxIter
        code=-1; % Too many iterations.
        break;
    end
    % Perform line search along p.
    [alpha,newX,newF]=linesearch(prob,veto,x,p,alphaMin,F,g'*p,mu,W,params);
	%alpha
	
    % Always update current point and residual. 
    x=newX;
    F=newF;

    % Store iteration trace and algorithm performance parameters.
    X=[X,x];
    alphas=[alphas alpha];

    if alpha==0
        code=-2; % Not sufficient descent.
        break;
    end
end


% Exit.

function [alpha,newx,newF]=linesearch(fun,vetoFun,x,p,alphaMin,F0,fp0,mu,W,params)
%LINESEARCH Perform Armijo linesearch with backtracking.
%
%[alpha,newx,newF]=linesearch(fun,vetoFun,x,p,alphaMin,F0,fp0,mu,W,params)
%fun      - name of the objective function.
%x        - current x value.
%p        - current search direction.
%alphaMin - shortest step-length allowed.
%F0       - residual function value at alpha=0.
%fp0      - gradient along line at alpha=0.
%mu       - Armijo constant.
%W        - weight matrix.
%params   - additional parameters to the objective function.
%alpha    - first step length that satisfied Armijos condition, or 0 if
%           no such alpha>=alphaMin is found.

% Calculate current objective function value.
f=0.5*F0'*W*F0;

% Start with full step length.
alpha=1.0;

% Continue while step length is long enough
while alpha>=alphaMin
    % Examine residual at proposed point.
    newx=x+alpha*p;
    newF=feval(fun,newx,params{:});
    newf=0.5*newF'*W*newF;

    % Is the reduction large enough?
    red=newf<f+mu*alpha*fp0;
    veto=false;
    if red && ~isempty(vetoFun)
        % Only calculate veto function if reduction is large enough.
        veto=feval(vetoFun,newx,params{:});
    end
	if red && ~veto
        % If reduction is large enough, we're done.
        return;
    else
        % Otherwise, try with half the step length.
        alpha=alpha/2;
    end
end

% No acceptable reduction found. Return alpha=0.
alpha=0;
newx=x;
newF=F0;
