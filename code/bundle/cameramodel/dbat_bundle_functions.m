%DBAT_BUNDLE_FUNCTIONS Overview of DBAT functions used in bundle computations.
%
%   Each DBAT function that is used in bundle computations has the
%   same basic structure. The typical use of each function is one
%   of the following (using a two-parameter function as an example):
%
%       1) V=NAME(P1, P2)
%       2) [V,dV]=NAME(P1, P2)
%       3) [V,dV]=NAME(P1, P2, C1, C2)
%
%   In the first case, the function value V is computed from the input
%   parameters. In the second case, the struct dV with the analytical
%   Jacobian with respect to each input parameter is also returned.
%   Each field is named after the corresponding formal parameter,
%   prepended by a 'd'. Thus, a parameter 'P' corresponds to the field
%   'dP'. In the third case, any supplied logical parameters that is
%   FALSE suppresses the computation of the corresponding Jacobian.
%   Any field corresponding to a non-computed Jacobian has en empty
%   value.
%
%   NOTE: All Jacobian with respect to a matrix-valued function F
%   and/or any matrix parameter A is understood to be with respect to
%   VEC(F(...)) and VEC(A), respectively.
%
%   During development and testing, two more, 'undocumented', uses are
%   likely:
%
%       4) [V,dV,dVn]=NAME(...)
%       5) NAME('SELFTEST',VERBOSE)
%
%   In the first case, the structure dVn will contain the numerical
%   approximations of the analytical Jacobians in dV. The numerical
%   approximation is performed by the function JACAPPROX. In the last
%   case, a self-test function is called that compares the analytical
%   and numerical Jacobians. If VERBOSE is supplied and TRUE, the
%   detailed result is written to the command window. Otherwise, the
%   function is quiet unless the computed differences are above a
%   hard-coded threshold. The latter case is useful for automated
%   testing of many functions.
%
%SEE ALSO: VEC, JACAPPROX.
