function [times,C,spLB,spLC]=time_icip(N,nIO,nEO,nOP,onlyDiag)
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
spLC=nnz(LC)/numel(LC);

extractClock=now;
extractTime=(extractClock-cholClock)*86400;

% Invert LA. Actually faster than LA\speye(size(LA)).
UA=inv(LA);

UAclock=now;
UAtime=(UAclock-extractClock)*86400;

LBUA=LB*UA; %#ok<MINV>

LBUAclock=now;
LBUAtime=(LBUAclock-UAclock)*86400;

% Computing LBUA directly may cause out of memory. Accept blocks of
% about 256MB.

blockSize=256*1024^2;

blockCols=floor(min(round(blockSize/8/size(LB,1)),size(LB,2))/3)*3

% Compute diagonal elements of U'*U. ud0 contains 3 elements per block
% [u11 u22 u33]
ud0=full(sum(UA.^2,1));

stopClock=now;
diagTime=(stopClock-LBUAclock)*86400;

lapClock=stopClock;

if onlyDiag
    offDiagTime=0;
    u12=[];
    u13=[];
    u23=[];
else
    % Extract staggered columns of U with stride 3
    ua1=UA(:,1:3:end);
    ua2=UA(:,2:3:end);
    ua3=UA(:,3:3:end);

    % Compute the (1,2), (1,3), (2,3) elements of each block. Each
    % vector will contain one element per block.
    u12=full(sum(ua1.*ua2,1));
    u13=full(sum(ua1.*ua3,1));
    u23=full(sum(ua2.*ua3,1));

    stopClock=now;
    offDiagTime=(stopClock-lapClock)*86400;

    lapClock=stopClock;
end

% Reset timers.
UBtime=0;

for base=1:blockCols:size(LBUA,2)
    % Columns in this block.
    ix=base:min(base+blockCols-1,size(UA,2));

    LBUAblk=LBUA(:,ix);
    UB=-LC\LBUAblk;

    stopClock=now;
    stopTime=(stopClock-lapClock)*86400;

    UBtime=UBtime+stopTime;
    lapClock=stopClock;
    
    % Each diagonal 3-by-3 block of
    % U'*U = [u11 u12 u13
    %         u12 u22 u23
    %         u13 u23 u33]
    
    % Update diagonal elements of U'*U.
    ud0(ix)=ud0(ix)+sum(UB.^2,1);

    stopClock=now;
    stopTime=(stopClock-lapClock)*86400;

    diagTime=diagTime+stopTime;
    lapClock=stopClock;

    if ~onlyDiag
        % Extract staggered columns of U with stride 3
        ub1=UB(:,1:3:end);
        ub2=UB(:,2:3:end);
        ub3=UB(:,3:3:end);

        ix3=(ix(1:3:end)-1)/3+1;
        % Compute the (1,2), (1,3), (2,3) elements of each block. Each
        % vector will contain one element per block.
        u12(ix3)=u12(ix3)+sum(ub1.*ub2,1);
        u13(ix3)=u13(ix3)+sum(ub1.*ub3,1);
        u23(ix3)=u23(ix3)+sum(ub2.*ub3,1);
        
        stopClock=now;
        stopTime=(stopClock-lapClock)*86400;

        offDiagTime=offDiagTime+stopTime;
        lapClock=stopClock;
    end
end

if onlyDiag
    C=spdiags(ud0(:),0,length(ud0),length(ud0));

    stopClock=now;
    combineTime=(stopClock-lapClock)*86400;
else
    % Create superdiagonal vectors
    udm1=reshape([u12;u23;zeros(size(u12))],[],1);
    udm2=reshape([u13;zeros(2,length(u13))],[],1);
    ud1=reshape([zeros(size(u12));u12;u23],[],1);
    ud2=reshape([zeros(2,length(u13));u13],[],1);

    % Preallocate place for all diagonals, including those
    % corresponding to unestimated OPs.

    ud=zeros(nOP,5);
    ud(:,1)=udm2;
    ud(:,2)=udm1;
    ud(:,3)=ud0;
    ud(:,4)=ud1;
    ud(:,5)=ud2;

    C=spdiags(ud,-2:2,size(ud,1),size(ud,1));

    stopClock=now;
    combineTime=(stopClock-lapClock)*86400;
end

stopClock=now;
totalTime=(stopClock-startClock)*86400;

times=[cholTime,extractTime,UAtime,LBUAtime,UBtime,diagTime,offDiagTime,combineTime,totalTime];
