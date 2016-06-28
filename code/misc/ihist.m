function h=ihist(v,size)
%IHIST Positive integer histogram.
%
%h=ihist(v[,size])
%v    - a data vector.
%size - size of histogram. Defaults to max(v).
%h    - a vector such that h(i) is the number of occurences of i in v.

% v1.0  96-01-24. Niclas Borlin, niclas@cs.umu.se.
% v1.01 96-11-22. Corrected a small bug that produced an error when
%                 non-integer data between 0 and 1 was supplied.
% v1.1  98-04-13. Now returns 0 if an empty vector is supplied.

if (nargin<2)
	size=max(v);
end

if (any(v<1))
	error('Input vector must be positive');
end

if (isempty(size))
	size=[1 1];
end
	
h=sparse(v,1,1,size,1);
