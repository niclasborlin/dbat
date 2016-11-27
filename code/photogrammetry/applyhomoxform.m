function y=applyhomoxform(A,x)
%APPLYHOMOXFORM Apply homogeneous transformation to non-homogeneous points.
%
%y=applyhomoxform(A,x)
%x - matrix with points as columns. Points may be homogeneous or non-homogeneous.
%A - homogeneous transformation matrix.
%y - points after transformation. y will be of same homogenity and size as x.


% Determine if x is homogeneous or not.
[m,n]=size(A);
isHomo=size(x,1)==m;
if ~isHomo
    % Add row of ones if not homogeneous.
    x=[x;ones(1,size(x,2))];
end

% Apply transformation.
y=A*x;

if ~isHomo
    % Normalize if input was not homogeneous.
    y=y(1:end-1,:)./repmat(y(end,:),m-1,1);
end
