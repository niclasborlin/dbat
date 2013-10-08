function y=applyhomoxform(A,x)
%APPLYHOMOXFORM Apply homogenous transformation to non-homogenous points.
%
%y=applyhomoxform(A,x)
%x - matrix with points as columns. Points may be homogenous or non-homogenous.
%A - homogenous transformation matrix.
%y - points after transformation. y will be of same homogenity and size as x.

% v1.0  2004-04-05. Niclas Borlin, niclas@cs.umu.se.

% Determine if x is homogenous or not.
[m,n]=size(A);
isHomo=size(x,1)==m;
if (~isHomo)
	% Add row of ones if not homogenous.
	x=[x;ones(1,size(x,2))];
end

% Apply transformation.
y=A*x;

if (~isHomo)
	% Normalize if input was not homogenous.
	y=y(1:end-1,:)./repmat(y(end,:),m-1,1);
end
