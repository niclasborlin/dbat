function CC=VectorizedCOPsparse(Lblocks,calc,useSparseLB)

useSparseLB=true;

% L = [LA, 0; LB, LC]. LA is square diagonal block corresponding to
% the OPs. LC is square diagonal block corresponding to the non-OPs.
% LB is rectangular subdiagonal block. LA is sparse block-diagonal
% with lower triangular 3-by-3 blocks. LB is dense. LC is dense lower
% triangular.
LA=Lblocks.LA;
LB=Lblocks.LB;
LC=Lblocks.LC;

% So, L*L'=J'*J.
% We want inv(J'*J) = inv(L*L') = inv(L')*inv(L) = U'*U.

% With U = [UA, 0; UB, UC] blocked as L, we are only interested in
% UA and UB.

if useSparseLB
    LB=Lblocks.LBsparse;
else
    LB=Lblocks.LB;
end

% Invert LA. Actually faster than LA\speye(size(LA)).
UA=inv(LA);

LBUA=LB*UA;
UB=-LC\LBUA;

if any(~calc)
    % Expand 
end

% Extract staggered columns of U with stride 3
ua1=UA(:,1:3:end);
ua2=UA(:,2:3:end);
ua3=UA(:,3:3:end);
ub1=UB(:,1:3:end);
ub2=UB(:,2:3:end);
ub3=UB(:,3:3:end);

% Each diagonal 3-by-3 block of
% U'*U = [u11 u12 u13
%         u12 u22 u23
%         u13 u23 u33]

% Pre-allocate vector for all elements, not just the one to be
% estimated.
ud0=zeros(size(calc));

% Compute diagonal elements of U'*U. ud contains 3 elements per
% block [u11 u22 u33]
ud0=(sum(UA.^2,1)+sum(UB.^2,1))';

% Compute the (1,2), (1,3), (2,3) elements of each block. Each
% vector will contain one element per block.
u12=sum(ua1.*ua2,1)+sum(ub1.*ub2,1);
u13=sum(ua1.*ua3,1)+sum(ub1.*ub3,1);
u23=sum(ua2.*ua3,1)+sum(ub2.*ub3,1);

% Create superdiagonal vectors
udm1=reshape([u12;u23;zeros(size(u12))],[],1);
udm2=reshape([u13;zeros(2,length(u13))],[],1);
ud1=reshape([zeros(size(u12));u12;u23],[],1);
ud2=reshape([zeros(2,length(u13));u13],[],1);

% Preallocate place for all diagonals, including those
% corresponding to unestimated OPs.

ud=zeros(numel(calc),5);
if ~isempty(udm2)
    ud(calc(:),1)=udm2;
    ud(calc(:),2)=udm1;
    ud(calc(:),3)=ud0;
    ud(calc(:),4)=ud1;
    ud(calc(:),5)=ud2;
end

CC=spdiags(ud,-2:2,size(ud,1),size(ud,1));
