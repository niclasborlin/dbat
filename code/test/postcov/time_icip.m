function [times,C,spLB]=time_icip(N,nIO,nEO,nOP,onlyDiag)
%Returned times are [chol,extract,UA,LBUA,UB,diag,offDiag,combine,total];

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

% Invert LA. Actually faster than LA\speye(size(LA)).
UA=inv(LA);

UAclock=now;
UAtime=(UAclock-extractClock)*86400;

LBUA=LB*UA; %#ok<MINV>

LBUAclock=now;
LBUAtime=(LBUAclock-UAclock)*86400;

UB=-LC\LBUA;

if 0
    blockCols=50000;
    blocks=ceil(size(LBUA,2)/blockCols);
    blockCols=ceil(size(LBUA,2)/blocks);

    tic
    UB2=zeros(size(LBUA));
    for k=1:blocks
        col1=(k-1)*blockCols+1;
        colN=min(col1+blockCols-1,size(LBUA,2));
        ix=col1:colN;
        Z=full(LBUA(:,ix));
        UB2(:,ix)=-LC\Z;
    end
    toc
end

UBclock=now;
UBtime=(UBclock-LBUAclock)*86400;

% Each diagonal 3-by-3 block of
% U'*U = [u11 u12 u13
%         u12 u22 u23
%         u13 u23 u33]

% Compute diagonal elements of U'*U. ud contains 3 elements per
% block [u11 u22 u33]
ud0=(full(sum(UA.^2,1))+sum(UB.^2,1))';

diagClock=now;
diagTime=(diagClock-UBclock)*86400;

if onlyDiag
    offDiagClock=diagClock;
    offDiagTime=(offDiagClock-diagClock)*86400;

    C=spdiags(ud0,0,length(ud0),length(ud0));

    combineClock=now;
    combineTime=(combineClock-offDiagClock)*86400;
else
    % Extract staggered columns of U with stride 3
    ua1=UA(:,1:3:end);
    ua2=UA(:,2:3:end);
    ua3=UA(:,3:3:end);
    ub1=UB(:,1:3:end);
    ub2=UB(:,2:3:end);
    ub3=UB(:,3:3:end);

    % Compute the (1,2), (1,3), (2,3) elements of each block. Each
    % vector will contain one element per block.
    u12=full(sum(ua1.*ua2,1))+sum(ub1.*ub2,1);
    u13=full(sum(ua1.*ua3,1))+sum(ub1.*ub3,1);
    u23=full(sum(ua2.*ua3,1))+sum(ub2.*ub3,1);

    % Create superdiagonal vectors
    udm1=reshape([u12;u23;zeros(size(u12))],[],1);
    udm2=reshape([u13;zeros(2,length(u13))],[],1);
    ud1=reshape([zeros(size(u12));u12;u23],[],1);
    ud2=reshape([zeros(2,length(u13));u13],[],1);

    % Preallocate place for all diagonals, including those
    % corresponding to unestimated OPs.

    ud=zeros(nOP,5);
    if ~isempty(udm2)
        ud(:,1)=udm2;
        ud(:,2)=udm1;
        ud(:,3)=ud0;
        ud(:,4)=ud1;
        ud(:,5)=ud2;
    end

    offDiagClock=now;
    offDiagTime=(offDiagClock-diagClock)*86400;

    C=spdiags(ud,-2:2,size(ud,1),size(ud,1));

    combineClock=now;
    combineTime=(combineClock-offDiagClock)*86400;
end

totalTime=(combineClock-startClock)*86400;

times=[cholTime,extractTime,UAtime,LBUAtime,UBtime,diagTime,offDiagTime,combineTime,totalTime];
