function v=abserr(A,B)
%ABSERR Absolute error, as measured by the Frobenius norm.
%
%   ABSERR(A,B), where A and B are same-sized arrays, returns the size
%   of the difference between A and B, as measured by the Frobenius
%   norm.

if any(size(A)~=size(B))
    error([mfilename,': argument size mismatch']);
end

v=norm(A-B,'fro');
