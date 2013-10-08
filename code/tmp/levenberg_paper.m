function [x,code,niter,X,deltas,rhos,F]=levenberg_paper(prob,vetoFun,x0,epsilon,maxIter,delta0,mu,eta,params)
%LEVENBERG Levenberg-Marquardts method for least squares problem.
%
%[x,code,niter,X,deltas,rhos]=levenberg(prob,vetoFun,x0,eps,maxIter,delta0,mu,eta,params)
%prob
%x0
%eps
%maxIter
%delta0
%mu
%eta
%x
%code
%iter
%X
%deltas
%rhos

x=x0;
X=x;
delta=delta0;
deltas=delta;
rhos=[];

niter=0;
[F,J]=feval(prob,x,params{:});
f=1/2*F'*F;
reduced=0;

while (1)
	[p,pGN]=dogleg(F,J,delta);
	%[p1,pGN1]=dogleglsq(F,full(J),delta);
    JpGN=J*pGN;
    Jp=J*p;
	if (norm(JpGN)<epsilon*norm(F) && all(p==pGN))
		break;
	end
	niter=niter+1;
	if (niter>maxIter)
		code=-1; % Too many iterations.
		return;
	end
	if (norm(p)==0)
		code=-2; % Not sufficient descent.
		return;
	end
	xNew=x+p;
	FNew=feval(prob,xNew,params{:});
	fNew=1/2*FNew'*FNew;
    if isempty(vetoFun)
        veto=false;
    else
        veto=feval(vetoFun,xNew,params{:});
    end
    if (veto)
        actual=-inf;
    else
        actual=f-fNew;
    end
	predicted=-F'*Jp-1/2*Jp'*Jp;
	rho=actual/predicted;
	rhos=[rhos rho];
	if rho<=mu
		reduced=reduced+1;
		if (reduced>maxIter)
			code=-1;
			break;
		end
		delta=delta/2;
        while delta>norm(pGN)
            % Keep adjusting delta until we're below norm(pGN). Otherwise
            % we would calculate the same steps over and over.
            delta=delta/2;
        end
		%niter=niter-1;
	else
		reduced=0;
		x=xNew;
        [F,J]=feval(prob,x,params{:});
		F=FNew;
		f=fNew;
		if rho>=eta
			delta=delta*2;
		end
		X=[X x];
		deltas=[deltas delta];
	end
end

code=0; % OK


function [p,pGN]=dogleg(F,J,delta)
%DOGLEGLSQ Perform a double dogleg step in the Levenberg-Marquardt method.
%
%[p,pGN]=dogleglsq(F,J,delta)
%F     - residual at current point.
%J     - Jacobian at current point.
%delta - current size of region of trust.
%p     - double dogleg search direction, |p|<=delta.
%pGN   - Gauss-Newton search direction.

% v1.0  99-11-09. Niclas Borlin, niclas@cs.umu.se.

debug=0;

% Calculate Gauss-Newton direction.
% Article formuation.
% pN=H\(-g)
% Gauss-Newton formulation.
% pGN=J\(-F)

% pGN=-inv(J'J)*J'*F
% Use singular value decomposition of J (will be reused below if G-N step
% too long).
H=J'*J;
g=J'*F;
pGN=H\(-g);

if (norm(pGN)<=delta)
	% Newton direction is short enough. Accept it.
	if (debug), disp('Gauss-Newton'); end
	p=pGN;
	return;
end

% Calculate the Cauchy Point.
% lambdaStar=g'*g/(g'*H*g);
% CP=-lambdaStar*g;

Hg=H*g;
lambdaStar=g'*g/(g'*Hg);
CP=-lambdaStar*g;

if (norm(CP)>delta)
	% Cauchy Point outside region. Use scaled gradient.
	if (debug), disp('CP'); end
	% p=-g/norm(g)*delta
	p=-g/norm(g)*delta;
	return;
end

% Calculate Nhat.
% gamma=norm(g)^4/((g'*H*g)*(g'*(H\g)));
gamma=(g'*g)^2/((g'*Hg)*(g'*(-pGN)));
eta=0.8*gamma+0.2;
Nhat=eta*pGN;

% Single dogleg.
Nhat=pGN;

if (norm(Nhat)<=delta)
	% Nhat is inside. Use scaled Newton.
	if (debug), disp('Scaled Gauss-Newton'); end
	p=pGN/norm(pGN)*delta;
	return;
end

% Find point on CP-nHat line with norm delta.
A=sum((CP-Nhat).^2);
B=sum(2*CP.*(Nhat-CP));
C=sum(CP.^2)-delta^2;

lambda=(-B+sqrt(B^2-4*A*C))/(2*A);

if (debug), disp('On line'); end
p=CP+lambda*(Nhat-CP);
