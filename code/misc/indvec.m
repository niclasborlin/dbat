function varargout=indvec(n,base)
%INDVEC Create index vectors.
%
%   I=INDVEC(N), where N>=0 is a scalar integer, returns the index vector
%   1:N.
%
%   [I,L]=INDVEC(N) also returns the last used index L=N.
%
%   [I1,I2,...,IM]=INDVEC(N), where N is an integer M-vector, returns the
%   vectors I1, I2, ..., IM, with consecutive indices I1=1:N(1),
%   I2=N(1)+(1:N(2)), etc.
%
%   [I1,I2,...,IM,L]=... returns L as the last used index L=sum(N).
%
%   ...=INDVEC(N,B) considers the indices 1:B to have been allocated, i.e.
%   returns the lowest index as B+1 instead of 1.
%
%   Example: Create consecutive index vectors of length 5, 4, 6, and 7,
%   respectively:
%
%       [i1,i2,i3,l]=indvec([5,4,6]);
%       i4=indvec(7,l);


if nargin<2, base=0; end

%varargout=cell(length(n)+1);
for i=1:length(n)
    varargout{i}=base+(1:n(i));
    base=base+n(i);
end

varargout{length(n)+1}=base;
