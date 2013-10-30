function [x,code,niter,X,lambdas,rhos,F]=levenberg_algebraic_paper(prob,vetoFun,x0,epsilon,maxIter,lambda0,minLambda,maxLambda,params)
%LEVENBERG Levenberg-Marquardts method for least squares problem.
%
%[x,code,niter,X,lambdas,rhos]=levenberg_algebraic_paper(prob,vetoFun,x0,epsilon,maxIter,lambda0,minLambda,maxLambda,params)
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
lambda=lambda0;
lambdas=lambda;
rhos=[];
prevLambda=lambda;

niter=0;
[F,J]=feval(prob,x,params{:});
f=1/2*F'*F;

while (1)
    JTJ=J'*J;
    JTF=J'*F;
    while niter<=maxIter
        % Solve for p.
        %log10(lambda)
        p=(JTJ+lambda*speye(size(J,2)))\(-JTF);
        niter=niter+1;
        Jp=J*p;
        xNew=x+p;
        FNew=feval(prob,xNew,params{:});
        fNew=1/2*FNew'*FNew;
        actual=f-fNew;
        predicted=-F'*Jp-1/2*Jp'*Jp;
        rho=actual/predicted;
        rhos=[rhos rho];
        if (isempty(vetoFun))
            veto=false;
        else
            veto=feval(vetoFun,xNew,params{:});
        end
        if (rho<0 || veto)
            % Bad step, try larger lambda.
            if (lambda==0)
                lambda=minLambda;
            else
                lambda=lambda*10;
            end
            lambdas=[lambdas lambda];
        else
            % Good step.
            x=xNew;
            [F,J]=feval(prob,x,params{:});
            f=fNew;
            lambda=lambda/10;
            if lambda<minLambda
                lambda=0;
            end
            X=[X x];
            lambdas=[lambdas lambda];
            break;
        end            
	end
        
	if (prevLambda==0 && norm(Jp)<epsilon*norm(F))
		break;
	end
    prevLambda=lambda;
	if (niter>maxIter)
		code=-1; % Too many iterations.
		return;
	end
end

code=0; % OK
