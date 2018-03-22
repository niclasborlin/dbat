function [x,code,n,final,T,rr,lambdas]=levenberg_marquardt(...
    resFun,vetoFun,x0,W,maxIter,termFun,doTrace,lambda0,lambdaMin,params)
%LEVENBERG_MARQUARDT Levenberg-Marquardt least squares adjustment algorithm.
%
%   [X,CODE,I]=LEVENBERG_MARQUARDT(RES,VETO,X0,N,TERM,TRACE,L0,MINL)
%   runs the Levenberg-Marquardt least squares adjustment algorithm
%   with weight matrix W on the problem with residual function RES and
%   with initial values X0. A maximum of N iterations are allowed. The
%   handle TERM should point to a termination function (see
%   GAUSS_NEWTON_ARMIJO) that will return TRUE if we are close enough
%   to the solution. The final estimate is returned in X. The damping
%   algorithm uses L0 as the initial lambda value. Any lambda value
%   below MINL is considered to be zero. In addition, if supplied and
%   non-empty, the VETO function is called to verify that the
%   suggested trial point is not invalid. The number of iteration I
%   and a success code (see GAUSS_NEWTON_ARMIJO) are also returned. If
%   TRACE is true, output sigma0 estimates at each iteration.
%
%   If the supplied L0 is negative, the initial lambda is calculated as
%   abs(L0)*trace(J0'*J0)/NN, where J0 is the Jacobian of the residual
%   function evaluated at X0, and NN is the number of unknowns. The same
%   applies for MINL.
%
%   [X,CODE,I,FINAL]=... also returns the struct FINAL with the final
%   step and the estimates of the weighted and unweighted residual
%   vector and Jacobian matrix. The final step are returned in the
%   field p. The weighted estimates are returned as fields weighted.r
%   and weighted.J, respectively, the unweighted as unweighted.r and
%   unweighted.J, respectively.
%
%   [X,CODE,I,FINAL,T,RR,LAMBDAS]=... returns the iteration trace as
%   successive columns in T, the successive estimates of sigma0 in RR and
%   the used damping values in LAMBDAS.
%
%   The function RES is assumed to return the residual function and its
%   jacobian when called [F,J]=feval(RES,X0).
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%     Nocedal, Wright (2006), "Numerical Optimization", 2nd ed.
%       Springer, Berlin, Germany. ISBN 978-0-387-40065-5.
%     Levenberg (1944), "A method for the solution of certain nonlinear
%       problems in least squares", Quarterly Journal of Applied
%       Mathematics, 2(2):164-168.
%     Marquardt (1963), "An algorithm for least squares estimation of
%       nonlinear parameters", SIAM Journal on Applied Mathematics,
%       11(2):431-441.
%
%See also: BUNDLE, GAUSS_MARKOV, GAUSS_NEWTON_ARMIJO,
%   LEVENBERG_MARQUARDT_POWELL.

% Initialize current estimate and iteration trace.
x=x0;

if nargout>4
    % Pre-allocate fixed block if trace is asked for.
    blockSize=50;
    T=nan(length(x),min(blockSize,maxIter+1));
end

% Iteration counter.
n=0;

% OK until signalled otherwise.
code=0;

% Compute Cholesky factor of weight matrix.
R=chol(W);

% Handle to weighted residual function.
wResFun=@(x)R*feval(resFun,x);

% Compute initial residual and Jacobian.
[s,K]=feval(resFun,x);
% Scale by Cholesky factor.
r=R*s;
J=R*K;
f=1/2*r'*r;
JTJ=J'*J;
JTr=J'*r;

% Residual norm trace.
rr=[];

% Compute real lambda0 if asked for.
if lambda0<0
    lambda0=abs(lambda0)*trace(JTJ)/size(J,2);
end

% Ditto for lambdaMin.
if lambdaMin<0
    lambdaMin=abs(lambdaMin)*trace(JTJ)/size(J,2);
end

% Initialize the damping parameter.
lambda=lambda0;

% Set to zero if below threshold.
if lambda<lambdaMin
    lambda=0;
end

% Store it.
lambdas=lambda;

% Damping from last successful step.
prevLambda=nan;

% Damping matrix.
I=speye(size(J,2));

while true
    % Stay in inner loop until a  better point is found or we run out of
    % iterations.
    while n<=maxIter
        % Solve for update p.
        p=(JTJ+lambda*I)\(-JTr);
        
        % Store current residual norm and used lambda value.
        rr(end+1)=sqrt(r'*r);

        % Verify the structural rank of the Jacobian. If it is less
        % than the number of columns, something is wrong.
        if n==0
            % A structural rank deficiency is a setup problems. Thus,
            % do the test only at the first iteration.
            srJ=sprank(J);
            if srJ<size(J,2)
                code=-4;
                p=nan(size(x));
                break;
            end
        end
        lambdas(end+1)=lambda;

        if doTrace
            if n==0
                fprintf(['Levenberg-Marquardt: iteration %d, ',...
                         'residual norm=%.2g\n'],n,rr(end));
            else
                fprintf(['Levenberg-Marquardt: iteration %d, ', ...
                         'residual norm=%.2g, lambda=%.2g\n'],n,rr(end),...
                        lambda);
            end
        end
    
        if nargout>4
            % Store iteration trace.
            if n+1>size(T,2)
                % Expand by blocksize if needed.
                T=[T,nan(length(x),blockSize)]; %#ok<AGROW>
            end
            T(:,n+1)=x;
        end

        % Update iteration count.
        n=n+1;
        
        % Pre-calculate J*p.
        Jp=J*p;

        % Compute weighted residual and objective function value at trial
        % point.
        t=x+p;
        rNew=feval(wResFun,t);
        fNew=1/2*rNew'*rNew;
        
        if fNew<f && ~isempty(vetoFun)
            % Call veto function only if we have a better point.
            fail=feval(vetoFun,t);
        else
            fail=false;
        end

        if fNew<f && ~fail
            % Good step, accept it.
            x=t;
            % Decrease damping in next iteration.
            lambda=lambda/10;
            % Switch to undamped if damping is small enough.
            if lambda<lambdaMin
                lambda=0;
            end

            % Evaluate unweighted residual and Jacobian at new point.
            [s,K]=feval(resFun,x);
            % Scale by Cholesky factor.
            r=R*s;
            J=R*K;
            f=1/2*r'*r;
            JTJ=J'*J;
            JTr=J'*r;
            
            % Leave inner loop since we found a better point.
            break;
        else
            % Bad step, discard t and increase damping.
            if lambda==0
                % Switch from undamped to minimum damping.
                lambda=lambdaMin;
            else
                lambda=lambda*10;
            end
        end
    end

    if code~=0
        break;
    end
    
    % We have either found a better point or run out of iterations.
    
    % Terminate with success if last step was without damping and we
    % satisfy the termination criteria.
    if prevLambda==0 && termFun(Jp,r)
        break;
    end
    
    % Remember lambda from last successful step.
    prevLambda=lambda;
    
    if n>maxIter
        code=-1; % Too many iterations.
        break;
    end
    
end

if nargout>3
    final=struct('unweighted',struct('r',s,'J',K),...
                 'weighted',struct('r',r,'J',J),...
                 'p',p);
end

if nargout>4
    % Store final point.
    T(:,n+1)=x;
end

rr(end+1)=sqrt(r'*r);

% Trim unused trace columns.
if nargout>4
    T=T(:,1:n+1);
end
