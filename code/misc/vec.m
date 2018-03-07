function v=vec(A)
%VEC Vectorization operation.
%
%   V=VEC(A) applies the vectorization operation on A, i.e. returns
%   the elements of A in column-major order.

v=A(:);
