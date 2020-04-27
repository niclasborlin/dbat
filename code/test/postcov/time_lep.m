function [times,C]=time_lep(N,nIO,nEO,nOP)
%Returned times are [total]
%
%References: Fraser 1987, "Limiting error propagation in network design".

if nIO>0
    nEO=nEO+nIO;
    nIO=0;
end

% N is partitioned EO-OP

% Macro to extract blocks of matrix. Blocksizes are M followed by N.
Blk11=@(A,m,n)A(1:m,1:m);
Blk12=@(A,m,n)A(1:m,m+1:m+n);
Blk22=@(A,m,n)A(m+1:m+n,m+1:m+n);

startClock=now;

D=Blk22(N,nEO,nOP);

% Invert D. Cannot use inv(D) or D\speye() since Matlab does not
% detect that D is block-diagonal.
DL=chol(D,'lower');
DLi=inv(DL);
C=DLi'*DLi;

stopClock=now;
totalTime=(stopClock-startClock)*86400;

times=[totalTime];
