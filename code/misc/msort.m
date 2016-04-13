function [y,i]=msort(A)
%MSORT Multiple level sort in ascending order.
%
%[y,i]=msort(A)
%y will be A sorted in ascending order. Ordering is determined by the
%  first column of A. Ties will be broken by looking at the second column.
%  Remaining ties will be broken by looking at the third column, etc.
%i will be an index vector such that y=A(i,:)

% v1.0  95-07-12. Niclas Borlin, niclas@cs.umu.se.

[m,n]=size(A);

% Start with no permutation.
i=1:m;

% Sort from right to left.
for j=n:-1:1
	% Sort this column.
	[ans,ix]=sort(A(i,j));
	% Modify permutation vector by result of this sort.
	i=i(ix);
end

y=A(i,:);
