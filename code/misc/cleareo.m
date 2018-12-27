function s=cleareo(s)
%CLEAREO Clear non-control object points in DBAT struct.
%
%   S=CLEAREO(S) clear all non-control external orientation
%   parameters, i.e. external orientation parameters that should be
%   estimated but do not have prior observations. The cleared
%   parameters are given a NaN value.

s.EO.val(s.bundle.est.EO & ~s.prior.EO.use)=NaN;
