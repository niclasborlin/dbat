function [x,code,n,final,T,rr,deltas,rhos,steps]=levenberg_marquardt_powell(...
    resFun,vetoFun,x0,W,maxIter,termFun,doTrace,delta0,mu,eta)
%LEVENBERG_MARQUARDT_POWELL Levenberg-Marquardt algorithm with Powell dogleg.
%
%   [X,CODE,I]=LEVENBERG_MARQUARDT_POWELL(RES,VETO,X0,W,N,TERM,TRACE,D0,MU,ETA)
%   runs the trust-region verions of the Levenberg-Marquardt least
%   squares adjustment algorithm with weight matrix W and Powell
%   dogleg on the problem with residual function RES and with initial
%   values X0. A maximum of N iterations are allowed. The handle TERM
%   should point to a termination function (see GAUSS_NEWTON_ARMIJO)
%   that will return TRUE if we are close enough to the solution. The
%   final estimate is returned in X. The damping algorithm uses D0 as
%   the initial delta value. The quality of each proposed step is
%   determined by the constants 0 < MU < ETA < 1. In addition, if
%   supplied and non-empty, the VETO function is called to verify that
%   the suggested trial point is not invalid. The number of iteration
%   I and a success code (see GAUSS_NEWTON_ARMIJO) are also returned.
%   If TRACE is true, output sigma0 estimates at each iteration.
%
%   [X,CODE,I,FINAL]=... also returns the struct FINAL with the final
%   step and the estimates of the weighted and unweighted residual
%   vector and Jacobian matrix. The final step are returned in the
%   field p. The weighted estimates are returned as fields weighted.r
%   and weighted.J, respectively, the unweighted as unweighted.r and
%   unweighted.J, respectively.
%
%   [X,CODE,I,FINAL,T,RR,DELTAS,RHOS,STEPS]=... returns the iteration trace as
%   successive columns in T, the successive estimates of sigma0 in RR, the
%   used damping values in DELTAS, the computed gain ratios in RHOS, and the
%   step types in STEPS. The step types are:
%     0 - Gauss-Newton,
%     2 - Cauchy (steepest descent),
%     1 - Interpolated between Gauss-Newton and Cauchy.
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
%       problems in least squares". Quarterly Journal of Applied
%       Mathematics, 2(2):164-168.
%     Marquardt (1963), "An algorithm for least squares estimation of
%       nonlinear parameters". SIAM Journal on Applied Mathematics,
%       11(2):431-441.
%     Powell (1970), "A Hybrid Method for Nonlinear Equations". In
%       "Numerical Methods for Nonlinear Algebraic Equations", (Ed.
%       Rabinowitz). Gordon and Breach Science, London UK:87-114.
%     Moré (1983), "Recent Developments in Algorithms and Software for
%       Trust Region Methods". In "Mathematical Programming - The State
%       of the Art" (Eds. Bachem, Grötschel, Korte), Springer, Berlin,
%       Germany: 258-287.
%
%See also: BUNDLE, GAUSS_MARKOV, GAUSS_NEWTON_ARMIJO, LEVENBERG_MARQUARDT.

% Initialize current estimate and iteration trace.
x=x0;

if nargout>4
    % Pre-allocate fixed block if trace is asked for.
    blockSize=50;
    T=nan(length(x),min(blockSize,maxIter+1));
    % Enter x0 as first column.
    T(:,1)=x0;
end

% Iteration counter.
n=0;

% OK until signalled otherwise.
code=0;

% Initialize the damping parameter.
delta=delta0;
deltas=[];

% Gain ratios.
rhos=[];

% Step types.
steps=[];

% Compute Cholesky factor of weight matrix.
R=chol(W);

% Handle to weighted residual function. Works for single-return
% call only. Used by linesearch.
wResFun=@(x)R*feval(resFun,x);

% Compute residual, Jacobian, and objective function value.
[s,K]=feval(resFun,x);
% Scale by Cholesky factor.
r=R*s;
J=R*K;
f=1/2*r'*r;

% Residual norm trace.
rr=[];

% Step type strings for trace output.
stepStr={'GN','IP','CP'};

while true
    % Store current residual norm.
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
    
    % Find search direction using the Powell single dogleg algorithm.
    [p,pGN,step]=dogleg(r,J,delta);

    % Store used lambda value.
    deltas(end+1)=delta;
    steps(end+1)=step;
    
    JpGN=J*pGN;
    Jp=J*p;

    if step==0 && termFun(J*pGN,r)
        % Call termination detection function. See BUNDLE for a discussion on
        % selection of termination detection functions. Only terminate
        % with success if last step was without damping, i.e. a full
        % G-N step.
        break;
    end

    % Evalutate residual and objective function value in trial point.
    t=x+p;
    rt=feval(wResFun,t);
    ft=1/2*rt'*rt;
    if isempty(vetoFun)
        veto=false;
    else
        veto=feval(vetoFun,t);
    end

    % Compare actual vs. predicted reduction.
    predicted=-r'*Jp-1/2*Jp'*Jp;
    actual=f-ft;

    % Gain ratio.
    rho=actual/predicted;
    rhos(end+1)=rho;

    if doTrace
        fprintf(['Levenberg-Marquardt-Powell: iteration %d, ',...
                 'residual norm=%.2g, delta=%.2g, step=%s, rho=%.1f\n'],n,...
                rr(end),delta,stepStr{step+1},rho);
    end
    
    if veto || rho<=mu
        % Point failed veto test or reduction was too poor.

        % Discard trial point, i.e. x is unchanged.
        
        % Reduce trust region size to half.
        delta=delta/2;
        
        % If necessary, shrink delta below norm(pGN). Otherwise we would
        % calculate and discard the same trial points multiple times.
        if delta>norm(pGN)
            % This replaces a while loop.
            delta=delta/pow2(ceil(log2(delta/norm(pGN))));
        end
    else
        % Accept new point.
        x=t;
        
        % Calculate residual, Jacobian, and objective function value at
        % new point.
        [s,K]=feval(resFun,x);
        % Scale by Cholesky factor.
        r=R*s;
        J=R*K;
        f=1/2*r'*r;

        % Increase trust region size if reduction was good.
        if rho>=eta
            delta=delta*2;
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
    
    if n>maxIter
        code=-1;
        break;
    end
end

if nargout>4
    % Store final point.
    T(:,n+1)=x;
end

if nargout>3
    final=struct('unweighted',struct('r',s,'J',K),...
                 'weighted',struct('r',r,'J',J),...
                 'p',p);
end

% Trim unused trace columns.
if nargout>4
    T=T(:,1:n);
end

function [p,pGN,step]=dogleg(r,J,delta)
%DOGLEGLSQ Perform a double dogleg step in the Levenberg-Marquardt method.
%
%[p,pGN,step]=dogleglsq(r,J,delta)
%r     - residual at current point.
%J     - Jacobian at current point.
%delta - current size of region of trust.
%p     - double dogleg search direction, |p|<=delta.
%pGN   - Gauss-Newton search direction.
%step  - step type, 0 - Gauss-Newton step,
%                   1 - Interpolated step.
%                   2 - Cauchy (Steepest descent) step,

% Calculate the Gauss-Newton direction.

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
Jn2=sum(J.^2,1);
Jn=sqrt(Jn2);
% Construct sparse diagonal scaling matrix.
D=sparse(1:length(Jn),1:length(Jn),1./Jn,length(Jn),length(Jn));
% Scale.
Js=J*D;
% Compute scaled Hessian (may be used later) and gradient.
Hs=Js'*Js;
gs=Js'*r;
% Solve scaled normal equations.
q=Hs\(-gs);
% Unscale solution.
pGN=D*q;

if norm(pGN)<=delta
    % Gauss-Newton direction is within region of trust. Accept it.
    p=pGN;
    step=0; % Signal GN step.
    return;
end

% Calculate the Cauchy Point.
%
%               g'*g
% lambdaStar = ------
%              g'*J'*J*g
%
% CP = -lambdaStar * g.
%
% Using scaled J and g: g=inv(D)*gs, J=Js*inv(D).
%
% g'*g = gs'*inv(D^2)*gs
%
% g'*J'*J*g = gs'*inv(D^2)*Js'*Js*inv(D^2)*gs=gs'*inv(D^2)*Hs*inv(D^2)*gs.
%
% CP = -lambdaStar * inv(D)*gs.

invD=sparse(1:length(Jn),1:length(Jn),Jn,length(Jn),length(Jn));
invD2=sparse(1:length(Jn),1:length(Jn),Jn2,length(Jn),length(Jn));
invD2gs=invD2*gs;
g=invD*gs;

lambdaStar=g'*g/(invD2gs'*Hs*invD2gs);

CP=-lambdaStar*g;

if norm(CP)>delta
    % Cauchy Point outside region of trust. Use scaled negative gradient.
    p=-g/norm(g)*delta;
    step=2; % Signal Cauchy step.
    return;
end

% Find intersection of line CP-pGN and circle with radius delta, i.e.
% find k such that norm(CP+k*(pGN-CP))=delta.

% Coefficients for second order equation.
A=sum((CP-pGN).^2);
B=sum(2*CP.*(pGN-CP));
C=sum(CP.^2)-delta^2;

% Solve for positive root.
k=(-B+sqrt(B^2-4*A*C))/(2*A);

% Point on circle.
p=CP+k*(pGN-CP);

% Signal interpolated step.
step=1;
