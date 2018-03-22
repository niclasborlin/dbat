function x=invblock(L,p,ix,method)
%INVBLOCK Compute one block of normal matrix inverse.
%
%   X=INVBLOCK(L,p,IX) computes the block M(IX,IX) of the inverse M of the
%   normal matrix N. L should have been computed as chol(N(p,p))', where p
%   is a permutation vector.


if isempty(ix)
    x=[];
    return;
end

if isnan(L(1,1))
    x=nan(nnz(ix));
end

% Create inverse permutation.
invP=zeros(size(p));
invP(p)=1:length(p);

switch method
  case 'verify'
    % Reconstruct N and solve via full inverse.
    LLT=L*L';
    N=LLT(invP,invP);
    invN=inv(N);
    x=invN(ix,ix);
  case 'sqrt'
    % Let E contain the wanted columns of I. Then our wanted block is
    % B=E'*(A\E)=E'*(inv(A)*E). With A=L*L',
    % B=E'*(inv(L*L')*E)=E'*inv(L')*inv(L)*E. With W=L\E=inv(L)*E, B=W'*W.

    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Run computation on multiple rows to enable profiling.
    v1=Ip(:,ix);
    v2=L\v1;
    x=v2'*v2;
  case 'sqrtsplit'
    % Let E contain the wanted columns of I. Then our wanted block is
    % B=E'*(A\E)=E'*(inv(A)*E). With A=L*L',
    % B=E'*(inv(L*L')*E)=E'*inv(L')*inv(L)*E. With W=L\E=inv(L)*E, B=W'*W.

    % Do a split of L into [L1,  0;
    %                       L21, L2]
    % L1 will generally be block-diagonal - use sparse storage.
    % L21 will generally be sparse - use sparse storage.
    % L2 will generally be dense - use full storage.
    
    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Find jump in row density.
    colDens=full(sum(L~=0,2))/size(L,1);
    [~,i]=max(diff(colDens));
    ix1=1:i;
    ix2=i+1:size(L,1);
    
    % Extract (1,1) block and (2,1) block as sparse matrices.
    L1=L(ix1,ix1);
    L21=L(ix2,ix1);
    % Extract (2,2) block as dense matrix.
    L2=full(L(ix2,ix2));

    % Permute the identity matrix.
    Ip=speye(size(L,2));
    Ip=Ip(p,:);

    % Run computation on multiple rows to enable profiling.
    % [L1    0  ] [ a ]   [ c ]
    % [L21   L2 ] [ b ] = [ d ]
    % a=L1\c
    % b=L2\(d-L21*a)
    cd=Ip(:,ix);
    c=cd(ix1,:);
    d=cd(ix2,:);
    a=L1\c;
    L21a=L21*a;
    dmL21a=d-L21a;
    b=L2\dmL21a;
    x=a'*a+b'*b;
  case 'direct'
    Ip=speye(size(L,2));
    Ip=Ip(p,:);
    A=L'\(L\Ip(:,ix));
    x=A(invP(ix),:);
  otherwise
    error('Bad method');
end

