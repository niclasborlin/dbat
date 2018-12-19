function [rms,res]=bundle_residuals(s,e)
%BUNDLE_RESIDUALS Point marking residuals after bundle run.
%
%   [RMS,RES]=BUNDLE_RESIDUALS(S,E) returns RMS as the pointwise residual
%   RMS in pixels for the whole project. The NUM_OP-by-NUM_IMAGES array RES
%   contains the pointwise RMS in pixels for each measured point that was
%   part of the bundle.
%
%See also: BUNDLE, BUNDLE_RESULT_FILE.

nPts=nnz(s.IP.vis);

% Compute residuals for each point.
ptRes=sqrt(sum(((1./s.IO.sensor.pxSize(:,s.IP.cam)).*...
                reshape(e.final.unweighted.r(1:nPts*2),2,[])).^2,1));

% Residual array has same sparsity as vis.
res=double(s.IP.vis);
res(s.IP.vis)=ptRes(s.IP.ix(s.IP.vis));

% RMS over all points.
rms=sqrt(mean(ptRes.^2));

