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
ptCols=full(s.IP.ix(s.IP.vis));

% Extract corresponding standard deviations.
stdIPpx=s.IP.std(:,ptCols);

% Standard deviations are given in pixels. The residuals are in mm, so
% scale the variance.
stdIPmm=stdIPpx.*s.IO.sensor.pxSize(:,s.IP.cam(ptCols));

% Convert to variance.
varIP=stdIPmm(:).^2;
% Variance of prior IO observations.
varIO=s.prior.IO.std(s.prior.IO.use).^2;
% Prior EO observations.
varEO=s.prior.EO.std(s.prior.EO.use).^2;
% Prior OP observations.
varOP=s.prior.OP.std(s.prior.OP.use).^2;

% Covariance matrix is currently diagonal. Pre-allocate diagonal and
% put variances in correct places.
d=nan(s.post.res.ix.n,1);
d(s.post.res.ix.IP)=varIP(:);
d(s.post.res.ix.IO)=varIO(:);
d(s.post.res.ix.EO)=varEO(:);
d(s.post.res.ix.OP)=varOP(:);

% Create covariance matrix.
Cobs=spdiags(d,0,s.post.res.ix.n,s.post.res.ix.n);

% Use the inverse as the weight matrix.
W=inv(Cobs);

