function [K,P,B,KC]=test_distortion_params(s,e)
%TEST_DISTORTION_PARAMS Chi-square test of lens and affine distortion parameters.
%
%   [K,P,B,KC]=test_distortion_params(S,E) computes p values from the
%   cumulative Chi-square distribution function for the lens and
%   affine distortion parameters. The structures S and E are given to
%   and returned by BUNDLE, respectively. The vector K returns the p
%   values for each radial distortion koefficient K1, K2, K3, treated
%   separately (DOF=1). The scalar P returns the p value for the P1,
%   P2 parameters together (DOF=2). The vector B returns the p values
%   for each affine distortion parameter B1, B2, treated separately
%   (DOF=1). NaN values are returned for unestimated parameters.
%
%See also: CUMCHI2, HIGH_IO_CORRELATIONS.


% Estimated IO values and their covariances.
x=s.IO;
CIO=bundle_cov(s,e,'CIO');
nCams=size(s.IO,2);
K=nan(s.nK,nCams);
KC=nan(s.nK,nCams);
P=nan(1,nCams);
B=nan(2,nCams);

if ~isempty(e.final.factorized) && e.final.factorized.fail
    % Tests are useless if factorization failed.
    return;
end

for j=find(s.IOunique)
    % Test radial coefficients individually.
    for i=1:s.nK
        % Chi-square statistic is (x-mu)'*inv(C)*(x-mu), where x is N(mu,C).
        ii=3+i;
        ix=sub2ind(size(x),ii,j);
        if s.estIO(ii,j)
            v=(x(ii,j)-0)^2/CIO(ix,ix);
            K(i,j)=cumchi2(v,1);
        end
    end
    % Also test radial coefficients together.
    for i=1:s.nK
        % Chi-square statistic is (x-mu)'*inv(C)*(x-mu), where x is N(mu,C).
        ii=3+(1:i);
        ix=sub2ind(size(x),ii,repmat(j,size(ii)));
        if all(s.estIO(ii,j))
            v=(x(ii,j)'-0)*(CIO(ix,ix)\(x(ii,j)-0));
            KC(i,j)=cumchi2(v,i);
        end
    end
end

% Test tangential coefficients together.
ii=(7:8)';
for j=find(s.IOunique)
    ix=sub2ind(size(x),ii,repmat(j,size(ii)));
    if all(s.estIO(ii,j))
        v=x(ii,j)'*(CIO(ix,ix)\x(ii,j));
        P(j)=cumchi2(v,2);
    end
end

% Test affine coefficients individually.
for j=find(s.IOunique)
    for i=1:size(B,1)
        % Chi-square statistic is (x-mu)'*inv(C)*(x-mu), where x is N(mu,C).
        ii=8+i;
        ix=sub2ind(size(x),ii,j);
        if s.estIO(ii,j)
            v=(x(ii,j)-0)^2/CIO(ix,ix);
            B(i,j)=cumchi2(v,1);
        end
    end
end
