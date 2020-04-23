function [times,C,spLB]=time_icip_mex(N,nIO,nEO,nOP,onlyDiag)
%Returned times are [chol,extract,mex,total];

%N is IO-EO-OP on entry.

p=[nIO+nEO+(1:nOP),nIO+(1:nEO),1:nIO];
N=N(p,p);

startClock=now;

L=chol(N,'lower');

cholClock=now;
cholTime=(cholClock-startClock)*86400;

% Extract blocks of L = [ A, 0; B, C].
% Diagonal OP block of L
LA=L(1:nOP,1:nOP);
% Diagonal non-OP block of L
LC=full(L(nOP+1:end,nOP+1:end));
% Subdiagonal block
LB=L(nOP+1:end,1:nOP);

spLB=nnz(LB)/numel(LB);

extractClock=now;
extractTime=(extractClock-cholClock)*86400;

if onlyDiag
    C=icpc_mex(LA,LB,LC,1);
else
    C=icpc_mex(LA,LB,LC,3);
end

mexClock=now;
mexTime=(mexClock-extractClock)*86400;

totalTime=(mexClock-startClock)*86400;

times=[cholTime,extractTime,mexTime,totalTime];
