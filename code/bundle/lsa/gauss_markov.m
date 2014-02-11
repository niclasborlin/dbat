function [x,code,n,f,J,T,rr]=gauss_markov(resFun,x0,maxIter,convTol,trace, ...
                                       sTest,params)
%GAUSS_MARKOV Gauss-Markov least squares adjustment algorithm.
%
%   [X,CODE,I]=GAUSS_MARKOV(RES,X0,N,TOL,TRACE,STEST,PARAMS) runs the
%   Gauss-Markov least squares adjustment algorithm on the problem with
%   residual function RES and with initial values X0. A maximum of N
%   iterations are allowed and the convergence tolerance is TOL. The final
%   estimate is returned in X. If STEST is true, the iterations are
%   terminated if a (near) singularity warning on the normal matrix is
%   issued. The number of iteration I and a success code (0 - OK, -1 - too
%   many iterations, -2 - matrix is singular) are also returned. If TRACE is
%   true, the sigma0 estimates are printed at each iteration.
%
%   [X,CODE,I,F,J]=... also returns the final estimates of the residual
%   vector F and Jacobian matrix J.
%
%   [X,CODE,I,F,J,T,RR]=... also returns the iteration trace as successive
%   columns in T and successive values of the residual norm in RR.
%
%   The function RES is assumed to return the residual function and its
%   Jacobian when called [F,J]=feval(RES,X0,PARAMS{:}), where the cell array
%   PARAMS contains any extra parameters for the residual function.
%
%   References:
%     Börlin, Grussenmeyer (2013), "Bundle Adjustment With and Without
%       Damping". Photogrammetric Record 28(144), pp. 396-415. DOI
%       10.1111/phor.12037.
%     McGlone et al. (2004), "Manual of Photogrammetry", 5th ed. American
%       Society for Photogrammetry and Remote Sensing, Bethesda, Maryland,
%       USA. ISBN 978-1570830716.
%
%See also: BUNDLE, GAUSS_NEWTON_ARMIJO, LEVENBERG_MARQUARDT,
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

while true
    % Calculate residual and Jacobian at current point.
    [f,J]=feval(resFun,x,params{:});

    rr(end+1)=sqrt(f'*f);
    if trace
        fprintf('Gauss-Markov: iteration %d, residual norm=%.2g\n',n,rr(end));
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

    % Terminate if angle between projected residual is smaller than
    % threshold. Warning! This test may be very strict on synthetic data
    % where norm(f) is close to zero at the solution.
    if norm(J*p)<=convTol*norm(f)
        % Converged.
        break
    end
    
    % Update iteration count.
    n=n+1;
    
    % Update estimate.
    x=x+p;

    if nargout>5
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

% Trim unused trace columns.
if nargout>5
    T=T(:,1:n+1);
end
