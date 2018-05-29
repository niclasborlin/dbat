function [x,code,n,final,T,rr,alphas]=gauss_newton_armijo(...
    resFun,vetoFun,x0,W,maxIter,termFun,trace,sTest,mu,alphaMin)
%GAUSS_NEWTON_ARMIJO Gauss-Newton-Armijo least squares adjustment algorithm.
%
%   [X,CODE,I]=GAUSS_NEWTON_ARMIJO(RES,VETO,X0,W,N,TERM,TRACE,STEST,MU,AMIN)
%   runs the Gauss-Newton-Armijo least squares adjustment algorithm
%   with weight matrix W on the problem with residual function RES and
%   with initial values X0. A maximum of N iterations are allowed. The
%   handle TERM should point to a termination function (see below)
%   that will return TRUE if we are close enough to the solution. The
%   final estimate is returned in X. If STEST is true, the iterations
%   are terminated if a (near) singularity warning on the normal
%   matrix is issued. The steplength algorithm uses the constant MU
%   for the Armijo condition and accepts steplengths down to AMIN. In
%   addition, if supplied and non-empty, the VETO function is called
%   to verify that the steplength does not generate an invalid point.
%   The number of iteration I and a success code (see below) are also
%   returned. If TRACE is true, sigma0 estimates are printed at each
%   iteration.
%
%   Termination function
%
%       The termination function should accept two vectors Jp and r
%       and return TRUE if they indicate that we are close enough to
%       the solution.
%
%       Standard termination functions are the relative termination
%       function
%
%             norm(Jp) <= tol*norm(r)
%
%       and the absolute termination function
%       
%             norm(r) <= tol,
%
%       where tol is a (small) constant.
%
%   Success codes
%
%       The following success codes are returned:
%
%        0  - OK
%       -1 - too many iterations,
%       -2 - normal matrix is singular
%       -3 - no alpha found by the line search
%       -4 - Jacobian is structurally rank deficient.
%
%   [X,CODE,I,FINAL]=... also returns the struct FINAL with the final
%   step and the estimates of the weighted and unweighted residual
%   vector and Jacobian matrix. The final step are returned in the
%   field p. The weighted estimates are returned as fields weighted.r
%   and weighted.J, respectively, the unweighted as unweighted.r and
%   unweighted.J, respectively.
%
%   [X,CODE,I,FINAL,T,RR,ALPHAS]=... returns the iteration trace as successive
%   columns in T, the successive estimates of sigma0 in RR, and the used
%   steplengths in ALPHAS.
%
%   The function RES is assumed to return the residual function and
%   its Jacobian when called [R,J]=feval(RES,X0).
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%     Nocedal, Wright (2006), "Numerical Optimization", 2nd ed.
%       Springer, Berlin, Germany. ISBN 978-0-387-40065-5.
%     Armijo (1966), "Minimization of Functions Having Lipschitz
%       Continuous First Partial Derivatives", Pacific Journal of
%       Mathematics, 16(1):1-3.
%
%See also: BUNDLE, GAUSS_MARKOV, LEVENBERG_MARQUARDT,
%   LEVENBERG_MARQUARDT_POWELL.

% Initialize current estimate and iteration trace.
x=x0;

if nargout>4
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

% Compute Cholesky factor of weight matrix.
R=chol(W);

% Handle to weighted residual function. Works for single-return
% call only. Used by linesearch.
wResFun=@(x)R*feval(resFun,x);

while true
    % Calculate unweighted residual and Jacobian at current point.
    [s,K]=feval(resFun,x);
    % [s,K,KK]=feval(resFun,x); abserr(K,KK), relerr(K,KK)
    % Scale by Cholesky factor.
    r=R*s;
    J=R*K;

    rr(end+1)=sqrt(r'*r); %#ok<AGROW>
    if trace
        if isempty(alphas)
            fprintf('Gauss-Newton-Armijo: iteration %d, residual norm=%.2g\n', ...
                    n,rr(end))
        else
            fprintf(['Gauss-Newton-Armijo: iteration %d, residual norm=%.2g, ' ...
                     'last alpha=%s\n'],n,rr(end),...
                    strtrim(rats(alphas(end))));
        end
    end

    % Verify the structural rank of the Jacobian. If it is less
    % than the number of columns, something is wrong.
    if n==0
        % A structural rank deficiency is a setup problems. Thus,
        % do the test only at the first iteration.
        srJ=sprank(J);
        if srJ<size(J,2)
            code=-4;
            p=nan(size(x));
            D=[]; Js=[]; % To avoid breaking 'final' structure.
            break;
        end
    end
    
    % Solve normal equations. Corresponds to p=(K'*W*K)\-(K'*W*s).
    
    % Do column scaling of the Jacobian to reduce the condition number.
    %
    % Original equation system:
    %
    %   J'*Jp=-J'*r.
    %
    % Apply scaling matrix D to the left of LHS, and RHS. Insert D*inv(D)
    % in the middle:
    %
    %   D*J'*J*D*inv(D)*p=-D*J'*r;
    %
    % Substitute q=inv(D)*p and solve for q in:
    % 
    %   D*J'*J*D*q=-D*J'*r;
    %
    % Finally, solve for p in
    %
    %   inv(D)*p=q, => p=D*q;

    % Column norms.
    Jn=sqrt(sum(J.^2,1));
    % Construct sparse diagonal scaling matrix.
    D=sparse(1:length(Jn),1:length(Jn),1./Jn,length(Jn),length(Jn));
    % Scale.
    Js=J*D;
    % Solve scaled normal equations.
    q=(Js'*Js)\-(Js'*r);
    % Unscale solution.
    p=D*q;

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
    
    % Call termination detection function. See BUNDLE for a
    % discussion on selection of termination detection functions.
    if termFun(Jp,r)
        % Converged.
        break
    end

    % Update iteration count.
    n=n+1;
    
    % Perform line search along p.
    [alpha,xNew,rNew]=linesearch(wResFun,vetoFun,x,p,alphaMin,r,r'*Jp,mu);
	
    % Always update current point and residual. 
    x=xNew;
    r=rNew;

    % Store iteration trace and algorithm performance parameters.
    alphas(end+1)=alpha; %#ok<AGROW>
    
    if nargout>4
        % Store iteration trace.
        if n+1>size(T,2)
            % Expand by blocksize if needed.
            T=[T,nan(length(x),blockSize)]; %#ok<AGROW>
        end
        T(:,n+1)=x;
    end

    if alpha==0
        % Not sufficient descent, stop.
        code=-3;
        % Ensure that length of residual norm vector matches length of alpha.
        rr(end+1)=rr(end); %#ok<AGROW>
        break;
    end
    
    if n>maxIter
        % Too many iterations, stop.
        code=-1; %
        % Ensure that length of residual norm vector matches length of  alpha.
        rr(end+1)=sqrt(r'*r); %#ok<AGROW>
        break;
    end

end

if nargout>3
    final=struct('unweighted',struct('r',s,'J',K),...
                 'weighted',struct('r',r,'J',J),...
                 'scaled',struct('D',D,'J',Js),...
                 'p',p);
end

% Trim unused trace columns.
if nargout>4
    T=T(:,1:n+1);
end


function [alpha,xNew,rNew]=linesearch(fun,veto,x,p,alphaMin,r0,fp0,mu)
%LINESEARCH Perform Armijo linesearch with backtracking.

% Calculate current objective function value.
f0=1/2*(r0'*r0);

% Always try full step length first.
alpha=1;

% Continue while step length is long enough
while alpha>=alphaMin
    % Examine residual at trial point.
    t=x+alpha*p;
    r=feval(fun,t);
    f=1/2*(r'*r);

    % Is the reduction large enough to satisfy the Armijo condition?
    redOK=f<f0+mu*alpha*fp0;

    if redOK && ~isempty(veto)
        % Call veto function only if reduction is large enough.
        fail=feval(veto,t);
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
rNew=r0;
