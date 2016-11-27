function [rms,res]=bundle_residuals(s,e)
%BUNDLE_RESIDUALS Point marking residuals after bundle run.
%
%   [RMS,RES]=BUNDLE_RESIDUALS(S,E) returns RMS as the pointwise residual
%   RMS in pixels for the whole project. The NUM_OP-by-NUM_IMAGES array RES
%   contains the pointwise RMS in pixels for each measured point that was
%   part of the bundle.
%
%See also: BUNDLE, BUNDLE_RESULT_FILE.


nPts=nnz(s.vis);

% Compute residuals for each point.
ptRes=sqrt(sum((diag(s.IO(end-1:end))*...
                reshape(e.final.unweighted.r(1:nPts*2),2,[])).^2,1));

% Residual array has same sparsity as vis.
res=double(s.vis);
res(s.vis)=ptRes(s.colPos(s.vis));

% RMS over all points.
rms=sqrt(mean(ptRes.^2));

