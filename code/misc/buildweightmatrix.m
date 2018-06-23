function W=buildweightmatrix(s)
%BUILDWEIGHTMATRIX Build weight matrix for the bundle.
%
%   W=BUILDWEIGHTMATRIX(S) computes the weight matrix W to be used
%   with the bundle from the DBAT structure S. The weight matrix is
%   the inverse of the observation covariance matrix. The observations
%   are all visible image observations and any prior IO, EO, or OP
%   observations.

% Image point variances. 

% Column indices into markPts that correspond to image measurements.
ptCols=full(s.colPos(s.vis));

% Extract corresponding standard deviations.
stdIPpx=s.markStd(:,ptCols);

% Standard deviations are given in pixels. The residuals are in mm, so
% scale the variance.
stdIPmm=stdIPpx./s.IO(end-1:end,s.ptCams(ptCols));

% Convert to variance.
varIP=stdIPmm(:).^2;
% Variance of prior IO observations.
varIO=s.prior.IOstd(s.useIOobs).^2;
% Prior EO observations.
varEO=s.prior.EOstd(s.useEOobs).^2;
% Prior OP observations.
varOP=s.prior.OPstd(s.useOPobs).^2;

% Covariance matrix is currently diagonal. Pre-allocate diagonal and
% put variances in correct places.
d=nan(s.residuals.ix.n,1);
d(s.residuals.ix.IP)=varIP(:);
d(s.residuals.ix.IO)=varIO(:);
d(s.residuals.ix.EO)=varEO(:);
d(s.residuals.ix.OP)=varOP(:);

% Create covariance matrix.
Cobs=spdiags(d,0,s.residuals.ix.n,s.residuals.ix.n);

% Use the inverse as the weight matrix.
W=inv(Cobs);

