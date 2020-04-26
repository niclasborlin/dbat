function [times,C]=time_classic_chol(N,nIO,nEO,nOP,sparseB,onlyDiag)
%Returned times are [prep, invert D, C1, C2, total]

% Reference: Wolf, Dewitt, Wilkinson 2013, "Elements of Photogrammetry
% with Application in GIS", Ch. 17-12.

if nIO>0
    times=nan(1,5);
    C=sparse(size(N,1),size(N,2));
    return;
end

% N is partitioned EO-OP

% Macro to extract blocks of matrix. Blocksizes are M followed by N.
Blk11=@(A,m,n)A(1:m,1:m);
Blk12=@(A,m,n)A(1:m,m+1:m+n);
Blk22=@(A,m,n)A(m+1:m+n,m+1:m+n);

startClock=now;

A=Blk11(N,nEO,nOP);
D=Blk22(N,nEO,nOP);
B=Blk12(N,nEO,nOP);

if sparseB
    B=sparse(B);
else
    B=full(B);
end

prepClock=now;
prepTime=(prepClock-startClock)*86400;

% Invert D. Cannot use inv(D) or D\speye() since Matlab does not
% detect that D is block-diagonal.
DL=chol(D,'lower');
DLi=inv(DL);
Di=DLi'*DLi;

diClock=now;
diTime=(diClock-prepClock)*86400;

% Compute C1, the EO post cov matrix
C1=inv(full(A-B*Di*B'));

c1Clock=now;
c1Time=(c1Clock-diClock)*86400;

% Compute C2, the OP post cov matrix
BDi=B*Di;

R1=chol(C1,'upper');

Z=R1*BDi;

if onlyDiag
    ZtZ=sum(Z.^2,1);
    CZ=spdiags(ZtZ(:),0,length(ZtZ),length(ZtZ));
    DiD=extractdiagblocks(Di,1);
    C2=DiD+CZ;
else
    error('x')
    C2f=diagblkouter(C1,BDi,3);
    C2=Di+C2f;
end

c2Clock=now;
c2Time=(c2Clock-c1Clock)*86400;

totalTime=(c2Clock-startClock)*86400;

allTimes=[prepTime,diTime,c1Time,c2Time,totalTime];

times=allTimes;
C=C2;
