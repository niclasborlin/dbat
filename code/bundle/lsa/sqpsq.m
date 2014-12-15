function [x,code,n,l,X,a,C,L,nus,r,J,A,rr]=sqpsq(fun,con,vetoFun,x0,maxIter,epsR,epsC,trace,singularTest,mu,nu0,alphaMin,params)
%SQPSQ Sequential quadratic programming for least-squares problems.
%
%[x,code,n,lambda,X,alphas,C,L,nu,r,J,A]=sqpsq(fun,con,x0,alphaMin,epsR,epsC,maxIter,params,nu0,mu)
%fun     - Name of residual/jacobian function.
%con     - Name of constraint/jacobian function.
%x0      - starting approximation.
%epsR    - convergence tolerance for ||Jp||/(1+||r||) (scalar).
%epsC    - convergence tolerance for |g| (scalar or vector).
%maxIter - maximum # of iterations.
%params  - cell array of parameters to pass along to fun and con.
%nu0     - starting penalty value for merit function.
%mu      - Armijo constant for line search.
%x       - The solution.
%n       - The needed number of iterations.
%code    - error code: 0 - OK, 1 - too many iterations.
%lambda  - vector of Lagrange multipliers at solution.
%X       - matrix with iterates as columns.
%alphas  - vector of used step lengths.
%C       - matrix with constrain function values at each iteration.
%L       - matrix with estimated Lagrange multipliers at each iteration.
%nu      - vector of used penalty weights.
%r       - residual vector at solution.
%J       - residual Jacobian at solution.
%A       - constraint Jacobian at solution.

n=0;
x=x0;
code=0;
X=[];
C=[];
L=[];
nu=nu0;
nus=nu;
a=[];
rr=[];

while (1)
    X=[X,x];
    % Evaluate functions and jacobians.
    [c,A]=feval(con,x,params{:});
    [r,J]=feval(fun,x,params{:});
    C=[C,c];
    rr(end+1)=sqrt(r'*r);

    % Construct the system matrix.
    Z=sparse(length(c),length(c));
    S=[-J'*J,A';A,Z];
    b=[J'*r;-c];
    % Solve the system matrix.
    pp=S\b;
    % Decompose solution vector.
    p=pp(1:length(x));
    l=pp(length(x)+[1:length(c)]);
    L=[L,l];

    Jp=J*p;
    if norm(Jp)<=epsR*(1+norm(r)) && all(abs(c)<=epsC)
        break;
    end

    % Calculate weight for merit function.
    oldNu=nu;
    curNu=selectweight(r,J,p,c,A,oldNu,0.1,0.1);
    nu=max(oldNu,curNu);

    nus=[nus,nu];

    % Calculate steplength.
    alpha=linesearch(fun,con,x,p,alphaMin,nu,r,c,-mu*((J'*r+nu*(A'*c))'*p),params);
    a=[a,alpha];
    if (alpha==0)
        code=2;
        break;
    end

    x=x+alpha*p;

    n=n+1;
    if n>maxIter
        code=1;
        break;
    end
end

function alpha=linesearch(fun,con,x,p,alphaMin,nu,oldR,oldC,reqRed,params)
%LINESEARCH Perform Armijo linesearch with backtracking on merit function.
%
%fun      - name of the residual function.
%con      - name of the constraint function.
%x        - current x value.
%p        - current search direction.
%alphaMin - shortest step-length allowed.
%nu       - penalty for merit function.
%oldR     - current residual.
%oldC     - current constraint residual.
%reqRed   - required reduction of the objective function.
%params   - additional parameters to the objective function.
%alpha    - first step length that satisfied Armijos condition, or 0 if
%           no such alpha>=alphaMin is found.

% Calculate current objective function value.
f=0.5*oldR'*oldR+0.5*nu*(oldC'*oldC);

% Start with full step length.
alpha=1.0;

% Continue while step length is long enough
while (alpha>=alphaMin)
    % Examine residual and constraint at proposed point.
    newx=x+alpha*p;
    newR=feval(fun,newx,params{:});
    newC=feval(con,newx,params{:});
    % Evaluate merit function.
    newf=0.5*newR'*newR+0.5*nu*(newC'*newC);
    
    % Calculate obtained reduction of function value.
    red=f-newf;
    
    if (red>=reqRed*alpha)
        % If reduction is large enough, we're done.
        return;
    else
        % Otherwise, try with half the step length.
        alpha=alpha/2;
    end
end

% No acceptable reduction found. Return alpha=0.
alpha=0;


function nu=selectweight(r,J,p,c,A,oldNu,delta1,delta2)
%SELECTWEIGHT Select constraint weight for the merit function.
%
%r      - residual vector.
%J      - jacobian.
%p      - search direction.
%c      - constraint vector.
%A      - constraint jacobian.
%oldNu  - weight used in last iteration.
%delta1 - lower tolerance on unit step length.
%delta2 - upper tolerance on unit step length.
%nu     - calculated weight such that
%         alphamin(nu)=1-delta1  if alphamin(oldNu)<1-delta1
%         alphamin(nu)=1+delta2  if alphamin(oldNu)>1+delta2
%         oldNu                  otherwise

%Calculate alphamin from old weight.

%phi(alpha)  =1/2*(nu||Ap alpha+c||^2+||L(Jp alpha+f)||^2)
%            =1/2*(nu(Ap alpha+c)^T(Ap alpha+c)+(Jp alpha+f)^TW(Jp alpha+f)
%phi'(alpha) =nu(c^TAp+alpha p^TA^TAp)+f^TWJp+alpha p^TJ^TWJp
%phi''(alpha)=nu p^TA^TAp+p^TJ^TWJp
%min value at -phi'(0)/phi''(0)

nu=oldNu;
Ap=A*p;
Jp=J*p;
cAp=c'*Ap;
rJp=r'*Jp;
Ap2=Ap'*Ap;
JpJp=Jp'*Jp;
	
phip0=nu*cAp+rJp;
phib0=nu*Ap2+JpJp;

% alpha=-(nu*cAp+rJp)/(nu*Ap2+JpJp)
alpha=-phip0/phib0;

if alpha<1-delta1
    % minimum alpha too small; adjust weight to make alphamin=1-delta1
    nu=(-rJp-(1-delta1)*JpJp)/(cAp+(1-delta1)*Ap2);
elseif alpha>1+delta2
    % minimum alpha too large; adjust weight to make alphamin=1+delta2
    nu=(-rJp-(1+delta1)*JpJp)/(cAp+(1+delta1)*Ap2);
end
