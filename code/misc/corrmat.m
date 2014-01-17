function [R,s]=corrmat(C)
%CORRMAT Calculate correlation matrix from covariance matrix.
%
%   CORRMAT(C), returns the correlation matrix corresponding to the
%   covariance matrix C.  C should be N-times-N and positive definite.
%
%   The output correlation matrix R contains correlation coefficients,
%   
%       R(i,j) = C(i,j)/sqrt(C(i,i)*C(j,j)).
%
%   The correlation coefficients will always be in the range [-1,1].
%
%   [R,s]=CORRMAT(C) also returns a vector of the standard deviations.
%
%See also COV, CORRCOEF, VAR, STD.

% $Id$

error(nargchk(nargin,1,1));

% Extract standard deviations.
s=sqrt(diag(C));
n=length(s);
% Construct sparse diagonal matrix with inverse.
Sinv=sparse(1:n,1:n,1./s,n,n);
% Calculate correlation coefficients.
R=Sinv*C*Sinv;

% Preserve any NaN's.
i=isnan(R);

% Force all elements to be between -1 and 1.
R=max(min(R,1),-1);
% Set diagonal elements to exactly 1.
R(speye(size(R))~=0)=1;

% Return any NaN's.
R(i)=NaN;
