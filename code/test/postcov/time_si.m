function [times,C]=time_si(N,nIO,nEO,nOP,swap)
%Returned times are [si, post, total]

%N is IO-EO-OP on entry.

if swap
    p=[nIO+nEO+(1:nOP),nIO+(1:nEO),1:nIO];
    N=N(p,p);
    % OP indices
    OPix=1:nOP;
else
    % OP indices
    OPix=nIO+nEO+(1:nOP);
end

startClock=now;

[Z,Zpattern]=sparseinv(N);

siClock=now;
siTime=(siClock-startClock)*86400;

C=extractdiagblocks(Z(OPix,OPix),3);

postClock=now;
postTime=(postClock-siClock)*86400;

totalTime=(postClock-startClock)*86400;

times=[siTime,postTime,totalTime];
