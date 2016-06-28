function [x,code,n,final,T,rr]=gauss_markov(resFun,x0,W,maxIter,convTol,trace,sTest)
%GAUSS_MARKOV Gauss-Markov least squares adjustment algorithm.
%
%   [X,CODE,I]=GAUSS_MARKOV(RES,X0,W,N,TOL,TRACE,STEST) runs the
%   Gauss-Markov least squares adjustment algorithm with weight matrix
%   W on the problem with residual function RES and with initial
%   values X0. A maximum of N iterations are allowed and the
%   convergence tolerance is TOL. The final estimate is returned in
%   X. If STEST is true, the iterations are terminated if a (near)
%   singularity warning on the normal matrix is issued. The number of
%   iteration I and a success code (0 - OK, -1 - too many iterations,
%   -2 - matrix is singular) are also returned. If TRACE is true, the
%   sigma0 estimates are printed at each iteration.
%
%   [X,CODE,I,FINAL]=... also returns the struct FINAL with the final
%   estimates of the weighted and unweighted residual vector and
%   Jacobian matrix. The weighted estimates are returned as fields
%   weighted.r and weighted.J, respectively, the unweighted as
%   unweighted.r and unweighted.J, respectively.
%
%   [X,CODE,I,FINAL,T,RR]=... also returns the iteration trace as successive
%   columns in T and successive values of the residual norm in RR.
%
%   The function RES is assumed to return the residual function and its
%   Jacobian when called [R,J]=feval(RES,X0).
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%     McGlone, Mikhail, Bethel, eds. (2004), "Manual of Photogrammetry",
%       5th ed., Chapter 11.1.3.4, pp. 786-788. American Society of
%       Photogrammetry and Remote Sensing.
%
%See also: BUNDLE, GAUSS_NEWTON_ARMIJO, LEVENBERG_MARQUARDT,
%   LEVENBERG_MARQUARDT_POWELL.

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

% Compute Cholesky factor of weight matrix.
R=chol(W);

while true
    % Calculate unweighted residual and Jacobian at current point.
    [s,K]=feval(resFun,x);
    % Scale by Cholesky factor.
    r=R*s;
    J=R*K;
    
    rr(end+1)=sqrt(r'*r);
    if trace
        fprintf('Gauss-Markov: iteration %d, residual norm=%.2g\n',n,rr(end));
    end
    
    % Solve normal equations. Corresponds to p=(K'*W*K)\-(K'*W*s).
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

    % Terminate if angle between projected residual is smaller than
    % threshold. Warning! This test may be very strict on synthetic data
    % where norm(r) is close to zero at the solution.
    if norm(J*p)<=convTol*norm(r)
        % Converged.
        break
    end
    
    % Update iteration count.
    n=n+1;
    
    % Update estimate.
    x=x+p;

    if nargout>4
        % Store iteration trace.
        if n+1>size(T,2)
            % Expand by blocksize if needed.
            T=[T,nan(length(x),blockSize)]; %#ok<AGROW>
        end
        T(:,n+1)=x;
    end

    % Terminate with error code if too many iterations.
    if n>maxIter
        code=-1; % Too many iterations.
        break;
    end
end

if nargout>3
    final=struct('unweighted',struct('r',s,'J',K),...
                 'weighted',struct('r',r,'J',J));
end

% Trim unused trace columns.
if nargout>4
    T=T(:,1:n+1);
end
