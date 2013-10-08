function [x,n,code,lambda,xx,alphas,cc,ll,nus,F,J]=sqp_paper(stub,veto,x0,nu0,mu,mi,epsf,epsg,W,params)
%SQP Sequential quadratic programming, damped version.
%
%[x,n,code,lambda,xx,alphas,cc,ll,nus]=sqp(stub,x0,nu0,mu,mi,epsf,epsg,W,params)
%stub    - stub name of problem. Residual/constraint functions should be
%          called stub_f and stub_c, respectively.
%x0      - starting approximation.
%nu0     - starting constraint weight.
%mi      - maximum # of iterations.
%epsf    - convergence tolerance for |Jp|/(1+|F|) (scalar).
%epsg    - convergence tolerance for |g| (scalar or vector).
%W       - weight matrix.
%params  - parameters to pass along to objective functions etc.
%x       - solution.
%n       - number of iterations.
%code    - error code: 0 - OK, 1 - too many iterations.
%lambda  - vector of Lagrange multipliers.
%res     - residual at solution.
%xx      - iteration trace of x.
%cc      - iteration trace of c.
%ll      - iteration trace of lambda.
%nus     - iteration trace of nu.

% v1.0  1999-11-05. Niclas Borlin, niclas@cs.umu.se.
% v1.1  2000-07-27. Now passes along parameters to functions.
% v1.01 2001-09-20. Made all matrices sparse.

fun=[stub '_f'];
con=[stub '_c'];

alphaMin=1e-3;
kappa=nu0*ones(4,1);

n=0;
x=x0;
code=0;
xx=[];
cc=[];
ll=[];
nu=nu0;
nus=nu;
alphas=[];

oldS=[];

M=inv(W);
R=chol(W);

while (1)
	xx=[xx,x];
	% Evaluate functions and jacobians.
	[c,A]=feval(con,x,params{:});
	[F,J]=feval(fun,x,params{:});
	cc=[cc,c];

	Z=sparse(length(c),length(c));
	S=[-J'*W*J,A';A,Z];
	b=[J'*W*F;-c];
	% Solve system matrix.
	pp=S\b;
	% Decompose solution vector.
	p=pp(1:length(x));
	lambda=pp(length(x)+[1:length(c)]);
	ll=[ll,lambda];

    Jp=J*p;
	if (abs(Jp'*F)<=norm(Jp)*norm(F)*epsf & all(abs(c)<=epsg))
		break;
	end

	if (nu0<0)
		% Undamped solution.
		alpha=1
		alphas=[alphas,alpha];
	else
		% Damped solution.
		
		% Calculate weight for merit function.
		oldNu=nu;
		curNu=selectweight(F,J,W,p,c,A,oldNu,0.1,0.1);
		% Cannot go below n highest weights.
        if (1)
            nu=max(curNu,kappa(1));
            if (nu~=oldNu)
                disp(sprintf(['Decreased nu=%g proposed. oldNu=%g, ',...
                              'kappa(1)=%g',curNu,oldNu,kappa(1)]));
            end
            % Update kappa queue if new weight is higher than old.
            if (nu>oldNu)
                kappa(1)=nu;
                kappa=sort(kappa);
            end
        else
            nu=max(oldNu,curNu);
            if (curNu<oldNu)
                nu=max(curNu,kappa(1));
                disp(sprintf(['Decreased nu=%g proposed. Selecting nu=%g ' ...
                              'instead'],curNu,nu));
                % Update kappa queue if new weight is higher than old.
                if (nu>oldNu)
                    kappa(1)=nu;
                    kappa=sort(kappa);
                end
            end
        end
		nus=[nus,nu];
		% Calculate steplength.
		alpha=linesearch(fun,con,veto,x,p,alphaMin,nu,F,c,-mu*(J'*(W*F)+nu*(A'*c))'*p,W,params);
		alphas=[alphas,alpha];
		% Take undamped step.
		if (alpha==0)
			code=2;
			break;
		end
	end
	x=x+alpha*p;
	if (any(imag(x)~=0) | ~all(isfinite(x)))
		code=3;
		break;
	end
	n=n+1;
	if (n>mi)
		code=1;
		break;
	end
end

function alpha=linesearch(fun,con,vetoFun,x,p,alphaMin,nu,oldF,oldC,reqRed,W,params)
%LINESEARCH Perform Armijo linesearch with backtracking.
%
%[alpha,newx,newF]=linesearch(fun,x,p,alphaMin,oldF,reqRed,W,params)
%fun      - name of the objective function
%x        - current x value
%p        - current search direction
%alphaMin - shortest step-length allowed
%reqRed   - required reduction of the objective function
%W        - weight matrix.
%params   - additional parameters to the objective function.
%alpha    - first step length that satisfied Armijos condition, or 0 if
%           no such alpha>=alphaMin is found.

% Calculate current objective function value.
f=0.5*oldF'*W*oldF+0.5*nu*(oldC'*oldC);

% Start with full step length.
alpha=1.0;

% Continue while step length is long enough
while (alpha>=alphaMin)
    % Examine residual at proposed point.
    newx=x+alpha*p;
    newF=feval(fun,newx,params{:});
	newC=feval(con,newx,params{:});
    newf=0.5*newF'*W*newF+0.5*nu*(newC'*newC);
    
    if (isempty(vetoFun))
        veto=false;
    else
        veto=feval(vetoFun,newx,params{:});
    end

    % Calculate obtained reduction of function value.
    red=f-newf;
    
    if (red>=reqRed*alpha & ~veto)
        % If reduction is large enough, we're done.
        return;
    else
        % Otherwise, try with half the step length.
        alpha=alpha/2;
    end
end

% No acceptable reduction found. Return alpha=0.
alpha=0;


function nu=selectweight(f,J,W,p,c,A,oldNu,delta1,delta2)
%SELECTWEIGHT Select constraint weight.
%
%nu=selectweight(f,J,W,p,c,A,oldNu,delta1,delta2)
%f      - residual vector.
%J      - jacobian.
%W      - weight matrix.
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

% v1.0  2003-06-19. Niclas Borlin, niclas@cs.umu.se.

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
fWJp=f'*(W*Jp);
Ap2=Ap'*Ap;
JpWJp=Jp'*W*Jp;
	
phip0=nu*cAp+fWJp;
phib0=nu*Ap2+JpWJp;

% alpha=-(nu*cAp+fWJp)/(nu*Ap2+JpWJp)
alpha=-phip0/phib0;

if (alpha<1-delta1)
	% minimum alpha too small; adjust weight to make alphamin=1-delta1
	nu=(-fWJp-(1-delta1)*JpWJp)/(cAp+(1-delta1)*Ap2);
elseif (alpha>1+delta2)
	% minimum alpha too large; adjust weight to make alphamin=1+delta2
	nu=(-fWJp-(1+delta1)*JpWJp)/(cAp+(1+delta1)*Ap2);
end
