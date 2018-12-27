function s=clearop(s)
%CLEAROP Clear non-control object points in DBAT struct.
%
%   S=CLEAROP(S) clear all non-control object point coordinates, i.e.
%   object points that should be estimated but do not have prior
%   observations. The cleared coordinates are given a NaN value.
%
%See also: SETCPT.

s.OP.val(s.bundle.est.OP & ~s.prior.OP.use)=NaN;
