function x=invblock(R,p,ix,method)
%INVBLOCK Compute one block of normal matrix inverse.
%
%   X=INVBLOCK(R,p,IX) computes the block M(IX,IX) of the inverse M of the
%   normal matrix N. R was computed as chol(N(p,p)), where p is a
%   permutation vector.

% $Id$

% Create inverse permutation.
invP=zeros(size(p));
invP(p)=1:length(p);

switch method
  case 'verify'
    % Reconstruct N and solve via full inverse.
    RTR=R'*R;
    N=RTR(invP,invP);
    invN=inv(N);
    x=invN(ix,ix);
  case 'split'
    % Do a split of R into [R1, R12;
    %                        0   R2]
    % R1 will generally be block-diagonal - use sparse storage.
    % R12, R2 will generally be dense - use full storage.
    
    % Permute the identity matrix.
    Ip=speye(size(R,2));
    Ip=Ip(p,:);

    % We are interested in the following rows of the solution.
    wIx=invP(ix);
    
    % Find jump in column density.
    colDens=full(sum(R~=0,1))/size(R,1);
    [~,i]=max(diff(colDens));
    ix1=1:i;
    ix2=i+1:size(R,1);
    
    % Which parts do our wanted rows belong to?
    wIx1=wIx<=i;
    wIx2=wIx>i;
    
    % Extract (1,1) block and (1,2) block as sparse matrices.
    R1=R(ix1,ix1);
    R12=R(ix1,ix2);
    % Extract (2,2) block as dense matrix.
    R2=full(R(ix2,ix2));

    % Run computation on multiple rows to enable profiling.
    v1=Ip(:,ix);
    v2=R'\v1;
    x2=R2\v2(ix2,:);
    
    if all(wIx2)
        % Wanted rows are all in x2, we are done.
        x=x2(wIx-i,:);
    else
        % Which is the first row we need?
        minR=min(wIx);
        
        % Only backsubstitute until the first row we need.
        v3=R12(minR:end,:)*x2;
        v4=v2(minR:i,:)-v3;
        x1=R1(minR:i,minR:i)\v4;
    
        xx=[zeros(minR-1,size(x1,2));x1;x2];
    
        x=xx(invP(ix),:);
    end
  case 'direct'
    Ip=speye(size(R,2));
    Ip=Ip(p,:);
    A=R\(R'\Ip(:,ix));
    x=A(invP(ix),:);
  otherwise
    error('Bad method');
end

