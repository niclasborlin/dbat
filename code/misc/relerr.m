function v=relerr(A,B)
%RELERR Relative error, as measured by the Frobenius norm.
%
%   RELERR(A,B), where A and B are same-sized arrays, returns the
%   relative size of the difference between A and B, as measured by
%   the Frobenius norm, normalized by the size of B.

if any(size(A)~=size(B))
    error([mfilename,': argument size mismatch']);
end

v=abserr(A,B)/norm(B,'fro');
