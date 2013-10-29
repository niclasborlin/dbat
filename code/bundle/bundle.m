function [s,X]=bundle(s)
%BUNDLE Run bundle adjustment iterations on a camera network.
%
%   [S,OK]=BUNDLE(S), where S is a struct returned by PROB2DBATSTRUCT, runs
%   the damped bundle adjustment on the camera network in the structure
%   S. The parameter values in S are used as initial values. The cIO, cEO,
%   cOP fields of S are used to indicate which parameters are free. OK is
%   returned as true if the bundle converged within the allowed number of
%   iterations. On return, the parameter values in S are updated if the
%   bundle converged.
%
%   ...=BUNDLE(S,DAMP) specifies which damping to use, 'none' (classic
%   bundle with no damping), 'GNA' (Gauss-Newton with Armijo linesearch),
%   'LM' (original Levenberg-Marquardt) , 'LMP' (Levenberg-Marquardt with
%   Powell dogleg).
%
%   [S,OK,S0]=... returns the sigma0 for the last iteration.
%
%   [S,OK,S0,X]=... returns the successive parameter estimates as columns
%   in X.
%
%   [S,OK,s0,X,CXX]=... returns the covariance matrix CXX of X, scaled by s0.
%
%
%References: Börlin, Grussenmeyer (2013), "Bundle Adjustment with and
%   without Damping", Photogrammetric Record 28(144), pp. XX-YY. DOI
%   10.1111/phor.12037.
%
%See also:

% $Id$
