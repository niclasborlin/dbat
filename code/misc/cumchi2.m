function y=cumchi2(x,n)
%CUMCHI2 Cumulative chi-square distribution function.
%
%   CUMCHI2(X,N) returns the cumulative chi-square distribution at the
%   values in the vector X with N degrees of freedom. N must be scalar.
%
%   Naive implementation that computes the cdf value via integration.


if ~isscalar(n)
    error('DBAT:CUMCHI2:invalidArgument','N must be scalar');
end
    
% Preallocate result vector.
y=zeros(size(x));

% Function to integrate.
fun=@(x) chi2(x,n);

% Slow loop version.
for i=1:length(x)
    if x(i)>n^2 && x(i)>100
        % Avoid lengthy integrations when unnecessary.
        y(i)=1;
    elseif x(i)>0
        % Integrate.
        y(i)=quad(fun,0,x(i));
    end
end

% Chi-square probability density function.
function y=chi2(x,n)

k=2^(-n/2)/gamma(n/2);

% Deal with negative X values separately to allow vectorized computation of
% rest.
neg=x<0;

if length(neg)
    y=zeros(size(x));
    x=x(~neg);
    y(~neg)=k*x.^(n/2-1).*exp(-x/2);
else
    y=k*x.^(n/2-1).*exp(-x/2);
end
