function [nrank,svals,indices,svals_err,stats] = spnrank(A,tol,nsvals,opts)
%SPNRANK   Numerical rank of dense or sparse matrices.
%   SPNRANK(A) provides an estimate of the number of linearly
%     independent rows or columns of a matrix A.
%   SPNRANK(A,TOL) is the number of singular values of A
%     that are larger than TOL.
%   SPNRANK(A) uses the default TOL = max(size(A)) * eps(norm of A), where
%     the norm of A is estimated using NORMEST_ERR, a minor modification of 
%     Matlab's NORMEST.
%   Note: SPNRANK returns -1 for NRANK if the  algorithm fails.
%
%   SPNRANK works for real or complex full matrices in Matlab 7.3 or higher
%      and for real (but not complex) sparse matrices in Matlab 7.5 or higher. 
%
%   [NRANK, SVALS] = SPNRANK(A, [], NSVALS), for the default tolerance, or
%   [NRANK, SVALS] = SPNRANK(A, TOL, NSVALS), for a user supplied TOL,
%   will return estimates of the NSVALS singular values of A closest to 
%   TOL that are larger than TOL and the NSVALS singular values of A
%   closest to TOL that are smaller than TOL.  If there are not NSVALS
%   singular values smaller or larger than TOL then as many as exist are
%   returned.  The entries in SVALS will be ordered largest to smallest.
%   If NSVALS is a vector with two components then SPNRANK returns NSVALS(1)
%   singular values greater than TOL and NSVALS(2) singular values less than
%   TOL, if this many exist.  If there is one output argument to SPNRANK
%   the default value of NSVALS is 0 and if there are two or more output
%   arguments to SPNRANK the default value of NSVALS is 3.
%
%   [NRANK, SVALS, INDICES] = SPNRANK(A, . . . ) also returns the indices
%   of the estimated singular values in SVALS.  SVALS(i) is an estimation
%   for singular value number INDICES(i) of A.
%
%   [NRANK, SVALS, INDICES, SVALS_ERR] = SPNRANK(A, [], NSVALS) or
%   [NRANK, SVALS, INDICES, SVALS_ERR] = SPNRANK(A, TOL, NSVALS) will
%   also return error bounds for the calculated singular values.
%   For each i, A (or, to be precise, a perturbation of A where the
%   perturbation is order of magnitude ||A|| eps where eps is relative
%   machine precision) is guaranteed to have a singular value in the
%   interval
%         [SVALS(i) - SVALS_ERR(i), SVALS(i) + SVALS_ERR(i)].
%   The bound does not guarantee that the true singular value of A
%   is the ith singular value but this is usually the case.
%
%   Note that if TOL lies in such an interval then the calculated numerical
%   rank, NRANK, may be incorrect (see STATS.flag_svals_err, below). This 
%   can occur, for example, if TOL is very close to a singular value of A.
%   To reduce SVALS_ERR one can decrease OPTS.tol_eigs, or increase 
%   OPTS.thresh (see below). Also changing TOL and NSVALS can affect 
%   SVALS_ERR.
%
%   SPNRANK(A,TOL,NSVALS,OPTS) specifies options:
%   OPTS.thresh - for sparse A OPTS.thresh is the threshold in the LDL 
%       calculation (see LDL) [ scalar in [0, 0.5] | {0.01} ]. Only
%       used in Matlab 7.6 or higher.
%   OPTS.tol_smax - the tolerance in estimating the norm of A (if
%       required)using NORMEST_ERR, a minor modification of Matlab's  
%       NORMEST [ scalar | {1.e-6} ]
%   The remaining fields are used only when NSVALS > 0.  In this case 
%       EIGS_MAXTIME is used to estimate the singular values of A near TOL. 
%       EIGS_MAXTIME is identical to Matlab's EIGS, except an option to 
%       limit the run time has been added.
%   OPTS.tol_eigs: stopping tolerance in EIGS  [scalar | {1.0e-6}]
%   OPTS.maxit: maximum number of iterations in EIGS [integer | {300}]
%   OPTS.disp: diagnostic information display level in EIGS and SPNRANK 
%              [{0} | 1 | 2]
%   OPTS.maxtime:  EIGS_MAXTIME is terminated after OPTS.maxtime seconds 
%      [scalar > 0 | {Inf (no time limit) } ]
%
%   [NRANK, SVALS, INDICES, SVALS_ERR, STATS] = SPNRANK(A, ...) will return
%   the structure STATS where
%   STATS.tol is the tolerance used to calculate the numerical rank
%   STATS.smax is an estimate of the two norm of A, if it is calculated
%   STATS.smax_err is a bound on | smax - (singular value  of A closest
%       to smax)|, if smax is calculated 
%   STATS.time_smax is the elapsed (wall clock) time to calculate smax
%   STATS.time_ldl is the elapsed time to run LDL
%   STATS.flag_ldl is 0 if LDL succeeded and 1 if LDL failed
%   STATS.message_ldl is the error message produced by LDL, if LDL fails
%   STATS.max_L is the largest magnitude element in L produced by LDL
%   STATS.time_eigs is the elapsed time to run EIGS_MAXTIME
%   STATS.flag_eigs is the flag returned by EIGS (if STATS.flag_eigs is 0
%       then all the eigenvalues converged; otherwise not all converged.)
%   STATS.time_svals_err is the time to calculate SVALS_ERR
%   STATS.flag_svals_err = 1 when  TOL lies within at least one interval
%       [SVALS(i) - SVALS_ERR(i), SVALS(i) + SVALS_ERR(i)] and the 
%       calculated NRANK may be incorrect.  STATS.flag_svals_err = 0 
%       is consistent with a correctly calculated numerical rank, NRANK.
%       STATS.flag_svals_err is NaN if any component of SVALS or SVALS_ERR
%       is not finite (which is possible when flag_eigs is 1).
%   Examples:
%       load west0479; 
%       A=west0479*west0479;
%       spnrank(A)
%       %or
%       [nrank, svals, indices, svals_err, stats]=spnrank(A);
%       %or
%       opts.tol_eigs = 1.e-10;
%       [nrank, svals, indices, svals_err, stats]=spnrank(A,[],[],opts);

%   The algorithm uses Sylvester's inertial theorem:
%      For a Hermitian matrix C if B = F * C * F'
%      for an invertible matrix F, then B and C have the same
%      number of eigenvalues greater than zero.
%   We apply this result to the Hermitian matrix
%   C = B - tol * I where B = [0 A ; A' 0 ] and use Matlab's
%   ldl decomposition: [L,D,P]=ldl(C). With this transformation 
%            C = P * L * D * L' * P' 
%   By Sylvester's theorem  C and D have the same number of positive
%   eigenvalues.  This is also the number of eigenvalues of B greater
%   than tol which is equal to the number of singular values of A 
%   greater than tol.
%
%   Mathematically, in exact arithmetic, SPNRANK will correctly
%   determine the rank of A.  In computer arithmetic SPNRANK
%   correctly determines the rank of the calculated LDL factorization,
%   which may have computer arithmetic errors.  However, in computer
%   arithmetic we have (see p. 218 of Accuracy and Stability of Numerical
%   Algorithms, 2nd ed. by Higham)
%               C + E = P * L * D * L' * P'
%   where P, L and D are the calculated factors and where
%     |E| <= p(m+n) ( |C| + P |L||D||L'| P' ) eps + O( eps^2).     (eqn 1)
%   Here eps is relative machine precision and p(x) is a linear
%   polynomial in x.  Since E is of magnitude proportional to eps,
%   in computer arithmetic SPNRANK will correctly determine the
%   rank of A except when tol is very close (in an interval of
%   magnitude proportional to eps ) to a singular value of A.
%
%   If TOL is O(norm(A)*eps) and A is singular then SPNRANK
%   may not determine the correct numerical rank since
%   TOL may be close to the zero singular values of A.
%   Therefore one should select TOL larger than some small
%   value.  For example if A is m by n:
%           TOL >= max(m,n)* eps * norm(A) 
%   appears to work well.  Note that Matlab's RANK also can be 
%   inaccurate of its tolerance is chosen too small.
%
%   As discussed earlier, SVALS and SVALS_ERR can be used to
%   provide an indication that TOL should be changed.
%
%   The error bounds in SVALS_ERR are based on Theorem 5.5, p. 205, 
%   of Applied Numerical Linear Algebra by James Demmel. Calculation
%   of the bounds directly uses C and does not use the calculated
%   LDL factorization.  Therefore the accuracy of the calculated
%   error bounds is not affected by the potential growth in errors
%   from the factor |L| |D| |L'| in (eqn 1) above.  However, as with
%   any computation in floating point arithmetic, potentially there
%   are errors in the calculations and, to be precise,  it is not
%   guaranteed that A has a singular value in
%       [SVALS(i) - SVALS_ERR(i), SVALS(i) + SVALS_ERR(i)] 
%   but it is guaranteed that A+E does where
%             ||E|| <= q(m,n) ||A|| eps                            (eqn 2)
%   and q(m,n) is a lower order polynomial in m and n.  If one 
%   characterizes the error in A using norms, (eqn 2) is the best one
%   can expect of a floating point computation involving matrices.

% Written by Leslie Foster 8/24/2008, modified 2/6/2009 Version 1.0
% Copyright 2008, Leslie Foster

%display_info = 1;     
display_info = 0;   %set to 1 to display start times of subroutines
% display_info is set to 1 below when opts.disp in input and opts.disp>=1

version_str = version;
idot = find(version_str == '.');
version_num = str2num( version_str(1:idot(1)+1) );

if (nargin == 1 || isempty(tol) )
    if ( nargin <= 3 || ~isfield(opts,'tol_smax') )
       tol_smax = 1.e-6;
    else
       tol_smax = opts.tol_smax;
    end
    time_smax = clock;
    [smax, cnt_smax, smax_err] = normest_err(A,tol_smax);
    time_smax = etime(clock,time_smax);
    % normest seems faster and safer than svds
    %    eg: [U,S,V,flag] = svds(A); smax = S(1,1);
    tol = max(size(A))*eps(smax);
    stats.tol = tol;
    stats.smax = smax; 
    stats.smax_err = smax_err;
    stats.time_smax = time_smax;
else
    stats.tol = tol;
    stats.smax = -1;
    stats.smax_err = -1;
    stats.time_smax = -1;
end
if ( nargin <= 2 || isempty(nsvals) )
   if (nargout >= 2)
      nsvals = 3;
   else
      nsvals = 0;
   end
end
if ( nargin <= 3 || ~isfield(opts,'thresh') )
   thresh = 0.01;  % the default value for ldl
else
   thresh = opts.thresh;
end
if ( max(nsvals) > 0 )
    if ( nargin <= 3 )
        bopts.tol = 1.0e-6;
        % Note that svals_err can be used to check that this
        % value is satisfactory.  The choice bopts.tol = eps,
        % the default value in EIGS, can have long run times.
        bopts.disp = 0;       
    else
        bopts = opts;
        if isfield( opts, 'tol_eigs')
            bopts = rmfield(bopts,'tol_eigs');
            bopts.tol = opts.tol_eigs;
        else
            bopts.tol = 1.0e-6;
        end
        if ~isfield(opts,'disp')
            bopts.disp = 0;
        end
    end
end
if nargin >= 4
    if isfield( opts, 'disp')
        if opts.disp >= 1
            display_info = 1;  % display info as calc proceeds
        end
    end
end

[m,n]= size(A);
minmn = min(m,n);

if  issparse(A) 
   C = [sparse(m,m) A;A' sparse(n,n)] - tol * speye(m+n,m+n);
else
   C = [zeros(m,m) A;A' zeros(n,n)] - tol * speye(m+n,m+n);
end

time_ldl = clock;
if ( display_info == 1 )
   disp(['time at start of ldl calc.: ', ...
       int2str(round(time_ldl(4:6))), ' hrs min sec'])
end

try
   if issparse(A)
       if ( version_num >= 7.6 )
          [L,D,p] = ldl(C,thresh,'lower','vector');           
       else
          [L,D,p] = ldl(C,'lower','vector');
          % this will throw an error, caught below, in Matlab 7.4 or earlier
       end
       L = (L.').' ;    % fix a problem with Matlab's sparse matrix storage
       D = (D.').' ;
   else
       [L,D,p] = ldl(C,'lower','vector');
       % this will throw an error, caught below, in Matlab 7.2 or earlier
   end
   time_ldl = etime(clock,time_ldl);
   stats.time_ldl = time_ldl;
   stats.flag_ldl = 0;
   max_L = full( max( max( abs( L(:) ) ) ) ) ;
   stats.max_L = max_L;
catch
   ME = lasterror;  
   stats.message_ldl = ME;
   stats.flag_ldl = 1;   
   time_ldl = etime(clock,time_ldl);
   stats.time_ldl = time_ldl;
   svals = [];            % since no svals calculated
   indices =[];
   svals_err = [];
   nrank =-1;             % flag if algorithm fails
   return
end

% D is block diagonal with 1 by 1 and 2 by 2 blocks on the diagonal.
% We can calculate the eigenvalues of D in a stable manner as follows:
d1 = full(diag(D));
d2 = full(diag(D,1));
i2 = find(d2~=0);
lam = d1;
cof = zeros(m+n,1);
cof(i2) = d1(i2) - d1(i2+1);
cof(i2)= 2* ((cof(i2)>=0) - (cof(i2)<0)) .* d2(i2) ./  ...
      ( abs(cof(i2))+ abs(complex( cof(i2), 2*d2(i2) ) ) );
% tricks here:
%    write to avoid potential subtractive cancellation and to construct
%    explicit orthogonal transformation
%       so to solve the eigenvalue problem for a 2 x 2 block [a b; b c] use: 
%       lam1 = a + b * 2 * sign(a - c) * b /(|a-c| + sqrt((a-c)^2 + 4 * b^2)
%       lam2 = c - b * 2 * sign(a - c) * b /(|a-c| + sqrt((a-c)^2 + 4 * b^2)
%    set sign using ((cof(i2)>=0) - (cof(i2)<0)) since Matlab's sign
%        has sign(0)=0 which is not correct here
%    use abs(complex( cof(i2), 2*d2(i2) ) )rather than 
%       sqrt(cof(i2).^2 + 4*d2(i2).^2) ) to avoid potential overflow
lam(i2) = d1(i2) + d2(i2) .* cof(i2);
lam(i2+1) = d1(i2+1) - d2(i2) .* cof(i2);
% construct orthogonal matrix V with D = V * (DV) * V' and DV diagonal
Vdiag = sqrt( 1.0 ./ (1 + cof .^ 2 ) );
Vdiag( i2 + 1 ) = Vdiag( i2 );
Vsuper = zeros(m+n,1);
Vsuper(i2) = - cof(i2) .* Vdiag(i2);
V = spdiags([-Vsuper,Vdiag,[0;Vsuper(1:end-1)]],-1:1,m+n,m+n) ;
DV = spdiags( lam, 0, m+n, m+n ) ;
nrank = length( find(lam>0) );

if ( sum(lam == 0) > 0 )
    warning('MATLAB:spnrank:TolNearExactSval',...
    ['[0 A ; A'' 0] - tol*I has a calculated eigenvalue equal 0\n' ...
            '         indicating that tol is near an exact singular ' ...
            ' value of A.\n         The algorithm may not converge unless' ...
            ' you try a new value for tol.\n']);
end
svals = [];               % for the case that no singular values are calculated
indices = [];
svals_err = [];
if ( max(nsvals) > 0 && nargout >= 2 )
   bopts.issym = 1;
   nC = length(C);
   if length(nsvals) == 1
       nsvals_l = min(nsvals,nrank);
       nsvals_s = min(nsvals,minmn - nrank);
   else
       nsvals_l = min(nsvals(1), nrank );
       nsvals_s = min(nsvals(2), minmn - nrank);
   end
   nsvals_tot = nsvals_l + nsvals_s;
   time_eigs = clock;
   if ( display_info == 1 )
       disp(['time at start of eigs: ', ...
           int2str(round(time_eigs(4:6))), ' hrs min sec'])
   end
   if ( nsvals_s == nsvals_l || nsvals_l == nsvals_s + 1 )
      [Veigs,Deigs,flag_eigs] = eigs_maxtime(@afunC, nC, nsvals_tot, ...
                  'be', bopts, L, V, DV);
   elseif ( nsvals_s == nsvals_l + 1 )
      [Veigs,Deigs,flag_eigs] = eigs_maxtime(@afunC, nC, nsvals_tot, ...
                  'be', bopts, L, V, -DV);
      Deigs = - Deigs;
   else
      flag_eigs = 0;
      Veigs = [];
      Deigs = [];
      if nsvals_l > 0
          [Veigs,Deigs,flag_eigs] = eigs_maxtime(@afunC, nC, nsvals_l, ...
                  'la', bopts, L, V, DV);
      end
      if nsvals_s > 0
         [Veigs1,Deigs1,flag1] = eigs_maxtime(@afunC, nC, nsvals_s, ...
                  'sa', bopts, L, V, DV);
         Veigs = [Veigs, Veigs1];
         Deigs = [ Deigs , sparse(size(Deigs,1),size(Deigs1,2)); ...
                   sparse(size(Deigs1,2),size(Deigs,1)), Deigs1 ];
         flag_eigs = max(flag_eigs,flag1);
      end
   end
   stats.time_eigs = etime(clock,time_eigs);
   stats.flag_eigs = flag_eigs;
   evals = diag(Deigs);
   svals = tol + 1.0 ./ evals;
   svals = max(svals,0*svals);    % assume that negative evals of B
                                  % are numerically 0
   [svals,isort] = sort(svals,'descend');
   indices = (nrank-nsvals_l+1):(nrank+nsvals_s);
   indices = indices';
   pt(p) = 1:length(p);
   Veigs = Veigs(pt,isort);
   
   if ( nargout >= 4 )
      time_svals_err = clock;
      if ( display_info == 1 )
         disp(['time at start of svals_err calc.: ', ...
             int2str(round(time_svals_err(4:6))), ' hrs min sec'])
      end
      svals_err= zeros(size(svals));
      for it = 1:length(svals)
         svals_err(it) = norm( C * Veigs(:,it) - (svals(it)-tol)*Veigs(:,it)) ;
         %bounds based on Theorem 5.5, p. 205, 
         %Applied Numerical Linear Algebra by James Demmel         
      end
      stats.time_svals_err = etime(clock,time_svals_err);
      if ( min(isfinite(svals)) == 1 && min(isfinite(svals_err)) == 1 )
         stats.flag_svals_err = max( abs(svals - tol ) <= svals_err ); 
      else
         stats.flag_svals_err = NaN;
      end
  end
end

end  %spnrank

%**************************************************************************
%-------------------------------------------------------------------------%
% NESTED FUNCTIONS:
%   (1)  afunC -- called by eigs_maxtime ( a minor variation of 
%        Matlab's eigs)
%   (2)  normest_err --  a small variation of Matlab's normest.  It 
%        returns an error bound on the estimated .
%   (3)  eigs_maxtime -- a variation of Matlab's eigs that allows
%        termination of the routine if it runs longer than a specified
%        time
%   (4)  Matlab's eigs contains many nested functions which are copied
%        here to make eigs_maxtime complete
%-------------------------------------------------------------------------%
%**************************************************************************

function y = afunC(x,L,V,DV)
% a function that evaluates S*x where S = inv(L*V*DV*V'*L') and
% L, V and DV are determined by SPNRANK. In SPNRANK the routine
% EIGS_MAXTIME (a minor variation of Matlab's EIGS) calls afunC
% in order to determine approximated singular values, when
% requested in SPNRANK.
% 
d=diag(DV);
y = L' \ ( V * (d .\ (V' * ( L \ x ) ) ) );

end %afunC

%**************************************************************************
%-------------------------------------------------------------------------%

function [e, cnt, nrm_err] = normest_err(S,tol)
%NORMEST_ERR Estimate the matrix 2-norm.
%   NORMEST_ERR(S) is an estimate of the 2-norm of the matrix S.
%   NORMEST_ERR(S,tol) uses relative error tol instead of 1.e-6.
%   [nrm,cnt] = NORMEST_ERR(..) also gives the number of iterations used.
%   [nrm, cnt, nrm_err] = NORMEST_ERR(..) returns nrm_err.  Here
%         | nrm - ( a singular value of A ) | <= nrm_err.
%   It is highly likely, but not guaranteed, that the singular
%   value of A is the largest singular value.
%
%   This function is intended primarily for sparse matrices,
%   although it works correctly and may be useful for large, full
%   matrices as well.  Use NORMEST_ERR when your problem is large
%   enough that NORM takes too long to compute and an approximate
%   norm is acceptable.
%
%   Class support for input S:
%      float: double, single
%
%   See also NORM, COND, RCOND, CONDEST.

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.14.4.2 $  $Date: 2005/05/31 16:31:10 $

%   Changes to estimate errors made by L. Foster 8/16/2008

if nargin < 2, tol = 1.e-6; end
x = sum(abs(S))';
cnt = 0;
normx = norm(x);
e = normx;
if e == 0, nrm_err = 0; return, end
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   x = x / normx;
   x0 = x;
   Sx = S*x;
   if nnz(Sx) == 0
      Sx = rand(size(Sx));
   end
   e = norm(Sx);
   x = S'*Sx;
   %w = norm(x - e^2 * x0);
   %disp([w e^2])
   normx = norm(x);
   cnt = cnt+1;
end
if nargout == 3
   %bounds based on Theorem 5.5, p. 205, 
   %Applied Numerical Linear Algebra by James Demmel    
   w = norm(x - e^2 * x0);
   nrm_err = w / e ;
   if ( w <= e^2 )
       nrm_err = w / ( e + sqrt( e^2 - w ) );
   end
end

end %normest_err

%**************************************************************************
%-------------------------------------------------------------------------%
function  varargout = eigs_maxtime(varargin)
%EIGS_MAXTIME is a minor variation of Matlab's EIGS that includes an option
%   to terminate the routine early if it runs longer than a specified
%   value (see OPTS.maxtime below).  For simplicity 'EIGS' rather than 
%   'EIGS_MAXTIME' is used in the remaining documentation.
%EIGS  Find a few eigenvalues and eigenvectors of a matrix using ARPACK
%   D = EIGS(A) returns a vector of A's 6 largest magnitude eigenvalues.
%   A must be square and should be large and sparse.
%
%   [V,D] = EIGS(A) returns a diagonal matrix D of A's 6 largest magnitude
%   eigenvalues and a matrix V whose columns are the corresponding
%   eigenvectors.
%
%   [V,D,FLAG] = EIGS(A) also returns a convergence flag. If FLAG is 0 then
%   all the eigenvalues converged; otherwise not all converged.
%
%   EIGS(A,B) solves the generalized eigenvalue problem A*V == B*V*D. B
%   must be symmetric (or Hermitian) positive definite and the same size as
%   A. EIGS(A,[],...) indicates the standard eigenvalue problem A*V == V*D.
%
%   EIGS(A,K) and EIGS(A,B,K) return the K largest magnitude eigenvalues.
%
%   EIGS(A,K,SIGMA) and EIGS(A,B,K,SIGMA) return K eigenvalues. If SIGMA is:
%      'LM' or 'SM' - Largest or Smallest Magnitude
%   For real symmetric problems, SIGMA may also be:
%      'LA' or 'SA' - Largest or Smallest Algebraic
%      'BE' - Both Ends, one more from high end if K is odd
%   For nonsymmetric and complex problems, SIGMA may also be:
%      'LR' or 'SR' - Largest or Smallest Real part
%      'LI' or 'SI' - Largest or Smallest Imaginary part
%   If SIGMA is a real or complex scalar including 0, EIGS finds the
%   eigenvalues closest to SIGMA. For scalar SIGMA, and when SIGMA = 'SM',
%   B need only be symmetric (or Hermitian) positive semi-definite since it
%   is not Cholesky factored as in the other cases.
%
%   EIGS(A,K,SIGMA,OPTS) and EIGS(A,B,K,SIGMA,OPTS) specify options:
%   OPTS.issym: symmetry of A or A-SIGMA*B represented by AFUN [{false} | true]
%   OPTS.isreal: complexity of A or A-SIGMA*B represented by AFUN [false | {true}]
%   OPTS.tol: convergence: Ritz estimate residual <= tol*NORM(A) [scalar | {eps}]
%   OPTS.maxit: maximum number of iterations [integer | {300}]
%   OPTS.p: number of Lanczos vectors: K+1<p<=N [integer | {2K}]
%   OPTS.v0: starting vector [N-by-1 vector | {randomly generated}]
%   OPTS.disp: diagnostic information display level [0 | {1} | 2]
%   OPTS.cholB: B is actually its Cholesky factor CHOL(B) [{false} | true]
%   OPTS.permB: sparse B is actually CHOL(B(permB,permB)) [permB | {1:N}]
%   OPTS.maxtime: terminate EIGS early (without completion) when run
%        time exceeds OPTS.maxtime [scalar | {Inf} ]
%   Use CHOL(B) instead of B when SIGMA is a string other than 'SM'.
%
%   EIGS(AFUN,N) accepts the function AFUN instead of the matrix A. AFUN is
%   a function handle and Y = AFUN(X) should return
%      A*X            if SIGMA is unspecified, or a string other than 'SM'
%      A\X            if SIGMA is 0 or 'SM'
%      (A-SIGMA*I)\X  if SIGMA is a nonzero scalar (standard problem)
%      (A-SIGMA*B)\X  if SIGMA is a nonzero scalar (generalized problem)
%   N is the size of A. The matrix A, A-SIGMA*I or A-SIGMA*B represented by
%   AFUN is assumed to be real and nonsymmetric unless specified otherwise
%   by OPTS.isreal and OPTS.issym. In all these EIGS syntaxes, EIGS(A,...)
%   may be replaced by EIGS(AFUN,N,...).
%
%   Example:
%      A = delsq(numgrid('C',15));  d1 = eigs(A,5,'SM');
%
%   Equivalently, if dnRk is the following one-line function:
%      %----------------------------%
%      function y = dnRk(x,R,k)
%      y = (delsq(numgrid(R,k))) \ x;
%      %----------------------------%
%
%      n = size(A,1);  opts.issym = 1;
%      d2 = eigs(@(x)dnRk(x,'C',15),n,5,'SM',opts);
%
%   See also EIG, SVDS, ARPACKC, FUNCTION_HANDLE.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.45.4.8 $  $Date: 2007/05/23 18:54:51 $

%   Changes to set maximum run time made by L. Foster 08/16/2008

%   EIGS provides the reverse communication interface to ARPACK library
%   routines. EIGS attempts to provide an interface for as many different
%   algorithms as possible. The reverse communication interfaces are
%   documented in the ARPACK Users' Guide, ISBN 0-89871-407-9.

cputms = zeros(5,1);
% code added by LVF:
t_start = clock; %start time for eigs
% end code added by LVF
t0 = cputime; % start timing pre-processing
      
% Process inputs and do error-checking
if (nargout > 3)
   error('MATLAB:eigs:TooManyOutputs', 'Too many output arguments.')
end

% Error check inputs and derive some information from them
[A,Amatrix,isrealprob,issymA,n,B,classAB,k,eigs_sigma,whch, ...
   sigma,tol,maxit,p,info,eigs_display,cholB,permB, ...
   resid,useeig,afunNargs,maxtime] = ...
   checkInputs(varargin{:});

% Now have enough information to do early return on cases EIGS does not
% handle. For these cases, use the full EIG code.
if useeig
   fullEig(nargout);
   return
end

if strcmp(eigs_sigma,'SM') || ~ischar(eigs_sigma)
   % eigs(A,B,k,scalarSigma) or eigs(A,B,k,'SM'), B may be []
   % Note: sigma must be real for [s,d]saupd and [s,d]naupd
   % If sigma is complex, even if A and B are both real, we use [c,z]naupd.
   % This means that mode=3 in [s,d]naupd, which has
   % OP = real(inv(A - sigma*M)*M) and B = M
   % reduces to the same OP as [s,d]saupd and [c,z]naupd.
   % A*x = lambda*M*x, M symmetric (positive) semi-definite
   % => OP = inv(A - sigma*M)*M and B = M
   % => shift-and-invert mode
   mode = 3;
elseif isempty(B)
   % eigs(A,k,stringSigma) or eigs(A,[],k,stringSigma), stringSigma~='SM'
   % A*x = lambda*x
   % => OP = A and B = I
   mode = 1;
else
   % eigs(A,B,k,stringSigma), stringSigma~='SM'
   % A*x = lambda*B*x
   % Since we can always Cholesky factor B, follow the advice of
   % Remark 3 in ARPACK Users' Guide, and do not use mode = 2.
   % Instead, use mode = 1 with OP(x) = R'\(A*(R\x)) and B = I
   % where R is B's upper triangular Cholesky factor: B = R'*R.
   % Finally, V = R\V returns the actual generalized eigenvectors of (A,B).
   mode = 1;
end

if cholB || ((mode == 1) && ~isempty(B))
   % The reordering permutation permB is [] unless B is sparse
   [RB,RBT,permB] = CHOLfactorB;
end

permAsB = [];
if (mode == 3) && Amatrix % need lu(A-sigma*B)
   % The reordering permutation permAsB is [] unless A-sigma*B is sparse
   [L,U,P,permAsB] = LUfactorAminusSigmaB;
end % if (mode == 3) && Amatrix

% Allocate outputs and ARPACK work variables
if isrealprob
   if issymA % real and symmetric
      if strcmp(classAB,'single')
         aupdfun = 'ssaupd';
         eupdfun = 'sseupd';
      else
         aupdfun = 'dsaupd';
         eupdfun = 'dseupd';
      end
      lworkl = int32(p*(p+8));
      d = zeros(k,1,classAB);
   else % real but not symmetric
      if strcmp(classAB,'single')
         aupdfun = 'snaupd';
         eupdfun = 'sneupd';
      else
         aupdfun = 'dnaupd';
         eupdfun = 'dneupd';
      end
      lworkl = int32(3*p*(p+2));
      workev = zeros(3*p,1,classAB);
      d = zeros(k+1,1,classAB);
      di = zeros(k+1,1,classAB);
   end
   v = zeros(n,p,classAB);
   workd = zeros(n,3,classAB);
   workl = zeros(lworkl,1,classAB);
else % complex
   if strcmp(classAB,'single')
      aupdfun = 'cnaupd';
      eupdfun = 'cneupd';
   else
      aupdfun = 'znaupd';
      eupdfun = 'zneupd';
   end
   zv = zeros(2*n*p,1,classAB);
   workd = complex(zeros(n,3,classAB));
   zworkd = zeros(2*numel(workd),1,classAB);
   lworkl = int32(3*p^2+5*p);
   workl = zeros(2*lworkl,1,classAB);
   workev = zeros(2*2*p,1,classAB);
   zd = zeros(2*(k+1),1,classAB);
   rwork = zeros(p,1,classAB);
end
ldv = int32(n);
ipntr = zeros(15,1,'int32');
ido = int32(0); % reverse communication parameter, initial value
if isempty(B) || (mode == 1)
   bmat = 'I'; % standard eigenvalue problem
else
   bmat = 'G'; % generalized eigenvalue problem
end
nev = int32(k); % number of eigenvalues requested
ncv = int32(p); % number of Lanczos vectors
iparam = zeros(11,1,'int32');
% iparam(1) = ishift = 1 ensures we are never asked to handle ido=3
iparam([1 3 7]) = [1 maxit mode];
select = zeros(p,1,'int32');

% To Do: Remove this error when ARPACKC supports singles
if strcmp(classAB,'single')
   error('MATLAB:eigs:single', ...
      'EIGS does not support single precision inputs.')
end

% The ARPACK routines return to EIGS many times per each iteration but we
% only want to display the Ritz values once per iteration (if opts.disp>0).
% Keep track of whether we've displayed this iteration yet in eigs_iter.
eigs_iter = 0;

cputms(1) = cputime - t0; % end timing pre-processing

% Iterate until ARPACK's reverse communication parameter ido says to stop
while (ido ~= 99)

   t0 = cputime; % start timing ARPACK calls **aupd

   if isrealprob
      arpackc( aupdfun, ido, ...
         bmat, int32(n), whch, nev, tol, resid, ncv, ...
         v, ldv, iparam, ipntr, workd, workl, lworkl, info );
   else
      % The FORTRAN ARPACK routine expects the complex input zworkd to have
      % real and imaginary parts interleaved, but the OP about to be
      % applied to workd expects it in MATLAB's complex representation with
      % separate real and imaginary parts. Thus we need both.
      zworkd(1:2:end-1) = real(workd);
      zworkd(2:2:end) = imag(workd);
      arpackc( aupdfun, ido, ...
         bmat, int32(n), whch, nev, tol, resid, ncv, ...
         zv, ldv, iparam, ipntr, zworkd, workl, lworkl, rwork, info );
      workd = reshape(complex(zworkd(1:2:end-1),zworkd(2:2:end)),[n,3]);
   end

   if (info < 0)
      error('MATLAB:eigs:ARPACKroutineError', ...
         'Error with ARPACK routine %s: info = %d', ...
         aupdfun,full(double(info)))
   end

   cputms(2) = cputms(2) + (cputime-t0); % end timing ARPACK calls **aupd
   t0 = cputime; % start timing MATLAB OP(X)

   % Compute which columns of workd ipntr references
   cols = checkIpntr;

   % The ARPACK reverse communication parameter ido tells EIGS what to do
   switch ido
      case {-1,1} % abs(ido)==1 => workd(:,col2) = OP*workd(:,col1)
         switch mode
            case 1 % mode==1 => OP(x) = K*x
               if isempty(B) % standard eigenvalue problem
                  % OP(x) = A*x
                  workd(:,cols(2)) = Amtimes(workd(:,cols(1)));
               else % generalized eigenvalue problem
                  % OP(x) = R'\(A*(R\x))
                  workd(:,cols(2)) = ...
                     RBTsolve(Amtimes(RBsolve(workd(:,cols(1)))));
               end
            case 3 % mode==3 => OP(x) = inv(A-sigma*B)*B*x
               if isempty(B) % standard eigenvalue problem
                  workd(:,cols(2)) = AminusSigmaBsolve(workd(:,cols(1)));
               else % generalized eigenvalue problem
                  switch ido
                     case -1
                        workd(:,cols(2)) = Bmtimes(workd(:,cols(1)));
                        workd(:,cols(2)) = ...
                           AminusSigmaBsolve(workd(:,cols(2)));
                     case 1
                        % mode==3 and ido==1:
                        % workd(:,col2) = inv(A-sigma*B)*B*x
                        % but B*x is already pre-computed in workd(:,col3)
                        workd(:,cols(2)) = ...
                           AminusSigmaBsolve(workd(:,cols(3)));
                     otherwise
                        error('MATLAB:eigs:UnknownRCP',...
                           'Unknown reverse communication parameter.')
                  end % switch ido (inner)
               end % if isempty(B)
            otherwise % mode is not 1 or 3
               error('MATLAB:eigs:UnknownMode','Unknown mode.')
         end % switch (mode)
      case 2 % ido==2 => workd(:,col2) = B*workd(:,col1)
         if (mode == 3)
            workd(:,cols(2)) = Bmtimes(workd(:,cols(1)));
         else
            error('MATLAB:eigs:UnknownMode','Unknown mode.')
         end
      case 3 % ido==3 => EIGS does not know how to compute shifts
         % setting iparam(1) = ishift = 1 ensures this never happens
         warning('MATLAB:eigs:WorklShiftsUnsupported', ...
            ['EIGS does not support computing the shifts in workl.' ...
            ' Returning immediately.'])
         ido = int32(99);
      case 99 % ido==99 => ARPACK is done
      otherwise
         error('MATLAB:eigs:UnknownReverseCommParamFromARPACK',...
            ['Unknown value of reverse communication parameter' ...
            ' returned from %s.'],aupdfun)

   end % switch ido (outer)

   cputms(3) = cputms(3) + (cputime-t0); % end timing MATLAB OP(X)

   if eigs_display
      displayRitzValues;
   end

   % code added by LVF
   if ( etime(clock,t_start) > maxtime )    
       ido = 99;
       if eigs_display
         disp(' ')
         disp(['elapsed time for eigs exceeded opts.maxtime of ', num2str(maxtime)])
         disp('  ****  exit eigs early  ****')
       end
   end
   % end code added by LVF
   
end % while (ido ~= 99)

t0 = cputime; % start timing post-processing

if (info < 0)
   error('MATLAB:eigs:ARPACKroutineError', ...
      'Error with ARPACK routine %s: info = %d',aupdfun,full(info));
end % if (info < 0)

if (nargout >= 2)
   rvec = int32(true); % compute eigenvectors
else
   rvec = int32(false); % do not compute eigenvectors
end

if isrealprob
   if issymA
      arpackc( eupdfun, rvec, 'A', select, ...
         d, v, ldv, sigma, ...
         bmat, int32(n), whch, nev, tol, resid, ncv, ...
         v, ldv, iparam, ipntr, workd, workl, lworkl, info );
      if strcmp(whch,'LM') || strcmp(whch,'LA')
         d = flipud(d);
         if (rvec == 1)
            v(:,1:k) = v(:,k:-1:1);
         end
      end
      if ((strcmp(whch,'SM') || strcmp(whch,'SA')) && (rvec == 0))
         d = flipud(d);
      end
   else
      % If sigma is complex, isrealprob=true and we use [c,z]neupd.
      % So use sigmar=sigma and sigmai=0 here in dneupd.
      arpackc( eupdfun, rvec, 'A', select, ...
         d, di, v, ldv, sigma, 0, workev, ...
         bmat, int32(n), whch, nev, tol, resid, ncv, ...
         v, ldv, iparam, ipntr, workd, workl, lworkl, info );
      d = complex(d,di);
      if rvec
         d(k+1) = [];
      else
         zind = find(d == 0);
         if isempty(zind)
            d = d(k+1:-1:2);
         else
            d(max(zind)) = [];
            d = flipud(d);
         end
      end
   end
else
   zsigma = [real(sigma); imag(sigma)];
   arpackc( eupdfun, rvec, 'A', select, ...
      zd, zv, ldv, zsigma, workev, ...
      bmat, int32(n), whch, nev, tol, resid, ncv, zv, ...
      ldv, iparam, ipntr, zworkd, workl, lworkl, ...
      rwork, info );
   if issymA
      d = zd(1:2:end-1);
   else
      d = complex(zd(1:2:end-1),zd(2:2:end));
   end
   v = reshape(complex(zv(1:2:end-1),zv(2:2:end)),[n p]);
end

flag = processEUPDinfo(nargin<3);

if (issymA) || (~isrealprob)
   if (nargout <= 1)
      if isrealprob
         varargout{1} = d;
      else
         varargout{1} = d(k:-1:1,1);
      end
   else
      varargout{1} = v(:,1:k);
      varargout{2} = diag(d(1:k,1));
      if (nargout >= 3)
         varargout{3} = flag;
      end
   end
else
   if (nargout <= 1)
      varargout{1} = d;
   else
      cplxd = find(di ~= 0);
      % complex conjugate pairs of eigenvalues occur together
      cplxd = cplxd(1:2:end);
      v(:,[cplxd cplxd+1]) = [complex(v(:,cplxd),v(:,cplxd+1)) ...
         complex(v(:,cplxd),-v(:,cplxd+1))];
      varargout{1} = v(:,1:k);
      varargout{2} = diag(d);
      if (nargout >= 3)
         varargout{3} = flag;
      end
   end
end

if (nargout >= 2) && (mode == 1) && ~isempty(B)
   varargout{1} = RBsolve(varargout{1});
end

cputms(4) = cputime-t0; % end timing post-processing

cputms(5) = sum(cputms(1:4)); % total time

if (eigs_display == 2)
   printTimings;
end

%-------------------------------------------------------------------------%
% Nested functions
%-------------------------------------------------------------------------%

% checkInputs error checks the inputs to EIGS and also derives some
%   variables from them:
% A may be a matrix or a function applying OP.
% Amatrix is true if A is a matrix, false if A is a function.
% isrealprob is true if all of A, B and sigma are real, false otherwise.
% issymA is true if A is symmetric, false otherwise.
% n is the size of (square) A and B.
% B is [] for the standard problem. Otherwise it may be one of B, CHOL(B)
%   or CHOL(B(permB,permB)).
% classAB is single if either A or B is single, otherwise double.
% k is the number of eigenvalues to be computed.
% eigs_sigma is the value for sigma passed in by the user, 'LM' if it was
%   unspecified. eigs_sigma may be either a string or a scalar value.
% whch is the ARPACK string corresponding to eigs_sigma and mode.
% sigma is the ARPACK scalar corresponding to eigs_sigma and mode.
% tol is the convergence tolerance.
% maxit is the maximum number of iterations.
% p is the number of Lanczos vectors.
% info is the start value, initialized to 1 or 0 to indicate whether to use
% resid as the start vector or not.
% eigs_display is true if Ritz values should be displayed, false otherwise.
% cholB is true if CHOL(B) was passed in instead of B, false otherwise.
% permB may be [], otherwise it is the permutation in CHOL(B(permB,permB)).
% resid is the start vector if specified and info=1, otherwise all zero.
% useeig is true if we need to use EIG instead of ARPACK, otherwise false.
% afunNargs is the range of EIGS' varargin that are to be passed as
%   trailing parameters to the function as in afun(X,P1,P2,...).
% maxtime is the maximum elapsed time allowed for eigs. (%LVF)
   function [A,Amatrix,isrealprob,issymA,n,B,classAB,k, ...
         eigs_sigma,whch,sigma,tol,maxit,p,info,eigs_display,cholB,...
         permB,resid,useeig,afunNargs,maxtime] = checkInputs(varargin)
      % Process inputs and do error-checking

      % Process the input A or the inputs AFUN and N
      % Start to derive some qualities (real, symmetric) about the problem
      if isfloat(varargin{1})
         A = varargin{1};
         Amatrix = true;
      else
         % By checking the function A with fcnchk, we can now use direct
         % function evaluation on the result, without resorting to feval
         A = fcnchk(varargin{1});
         Amatrix = false;
      end
      % isrealprob = isreal(A) && isreal(B) && isreal(sigma)
      isrealprob = true;
      issymA = false;
      if Amatrix
         isrealprob = isreal(A);
         issymA = ishermitian(A);
         [m,n] = size(A);
         if (m ~= n)
            error('MATLAB:eigs:NonSquareMatrixOrFunction',...
               'A must be a square matrix or a function.')
         end
      else
         n = varargin{2};
         nstr = 'Size of problem, ''n'', must be a positive integer.';
         if ~isscalar(n) || ~isreal(n)
            error('MATLAB:eigs:NonPosIntSize', nstr)
         end
         if issparse(n)
            n = full(n);
         end
         if (round(n) ~= n)
            warning('MATLAB:eigs:NonPosIntSize',['%s\n         ' ...
               'Rounding input size.'],nstr)
            n = round(n);
         end
      end

      % Process the input B and derive the class of the problem.
      % Is B present in the eigs call or not?
      Bpresent = true;
      Bstr = ['Generalized matrix B must be the same size as A and' ...
         ' either a symmetric positive (semi-)definite matrix or' ...
         ' its Cholesky factor.'];
      if (nargin < (3-Amatrix))
         B = [];
         Bpresent = false;
      else
         % Is the next input B or K?
         B = varargin{3-Amatrix};
         if ~isempty(B) % allow eigs(A,[],k,sigma,opts);
            if isscalar(B)
               if n ~= 1
                  % this input is really K and B is not specified
                  B = [];
                  Bpresent = false;
               else
                  % This input could be B or K.
                  % If A is scalar, then the only valid value for k is 1.
                  % So if this input is scalar, let it be B, namely
                  % eigs(4,2,...) assumes A=4, B=2, NOT A=4, k=2
                  if ~isnumeric(B)
                     error('MATLAB:eigs:BsizeMismatchAorNotSPDorNotChol', Bstr);
                  end
                  % Unless, of course, the scalar is 1, in which case 
                  % assume the that it is meant to be K.
                  if (B == 1) && ((Amatrix && nargin <= 3) || ...
                         (~Amatrix && nargin <= 4))
                      B = [];
                      Bpresent = false;
                  elseif ~isfloat(B)
                     error('MATLAB:eigs:BsizeMismatchAorNotSPDorNotChol', Bstr);
                  end
               end
            else
               % B is a not a scalar.
               if ~isfloat(B) || ~isequal(size(B),[n,n])
                  error('MATLAB:eigs:BsizeMismatchAorNotSPDorNotChol', Bstr);
               end
               isrealprob = isrealprob && isreal(B);
            end
         end
      end
      % ARPACK can only handle homogeneous inputs
      if Amatrix
         classAB = superiorfloat(A,B);
         A = cast(A,classAB);
         B = cast(B,classAB);
      else
         if ~isempty(B)
            classAB = class(B);
         else
            classAB = 'double';
         end
      end
      
      % argOffset tells us where to get the eigs inputs K, SIGMA and OPTS.
      % If A is really the function afun, then it also helps us find the
      % trailing parameters in eigs(afun,n,[B],k,sigma,opts,P1,P2,...)
      % Values of argOffset:
      %  0: Amatrix is false and Bpresent is true:
      %     eigs(afun,n,B,k,sigma,opts,P1,P2,...)
      %  1: Amatrix and Bpresent are both true, or both false
      %     eigs(A,B,k,sigma,opts)
      %     eigs(afun,n,k,sigma,opts,P1,P2,...)
      %  2: Amatrix is true and Bpresent is false:
      %     eigs(A,k,sigma,opts)
      argOffset = Amatrix + ~Bpresent;

      if Amatrix && ((nargin - Bpresent)>4)
         error('MATLAB:eigs:TooManyInputs', 'Too many inputs.')
      end

      % Process the input K.
      if (nargin < (4-argOffset))
         k = min(n,6);
      else
         k = varargin{4-argOffset};
         kstr = ['Number of eigenvalues requested, k, must be a' ...
            ' positive integer <= n.'];
         if ~isnumeric(k) || ~isscalar(k) || ~isreal(k) || (k>n)
            error('MATLAB:eigs:NonIntegerEigQty', kstr)
         end
         if issparse(k)
            k = full(k);
         end
         if (round(k) ~= k)
            warning('MATLAB:eigs:NonIntegerEigQty',['%s\n         ' ...
               'Rounding number of eigenvalues.'],kstr)
            k = round(k);
         end
      end

      % Process the input SIGMA and derive ARPACK values whch and sigma.
      % eigs_sigma is the value documented in the help as "SIGMA" that is
      % passed in to EIGS. eigs_sigma may be either a scalar, including 0,
      % or a string, including 'SM'.
      % In ARPACK, eigs_sigma corresponds to two variables:
      % 1.  which, called "whch" to avoid conflict with MATLAB's function
      % 2.  sigma
      % whch is always a string. sigma is always a scalar.
      % Valid combinations are shown below. Note eigs_sigma = 0/'SM' has
      % the same sigma/whch values as eigs_sigma='LM' (default) so these
      % must be distinguished by the mode.
      % eigs_sigma = 'SM' or 0 => sigma = 0, whch = 'LM' (mode=3)
      % eigs_sigma is a string not 'SM' => sigma = 0, whch = eigs_sigma (mode=1)
      % eigs_sigma is any scalar => sigma = eigs_sigma, whch = 'LM'
      % (mode=1)
      whchstr = 'Eigenvalue range sigma must be a valid 2-element string.';
      if (nargin < (5-argOffset))
         % default: eigs 'LM' => ARPACK which='LM', sigma=0
         eigs_sigma = 'LM';
         whch = 'LM';
         sigma = 0;
      else
         eigs_sigma = varargin{5-argOffset};
         if ischar(eigs_sigma)
            % eigs(string) => ARPACK which=string, sigma=0
            if ~isequal(size(eigs_sigma),[1,2])
               error('MATLAB:eigs:EigenvalueRangeNotValid', ...
                  [whchstr '\nFor real symmetric A, the' ...
                  ' choices are ''%s'', ''%s'', ''%s'', ''%s'' or ''%s''.' ...
                  '\nFor non-symmetric or complex' ...
                  ' A, the choices are ''%s'', ''%s'', ''%s'', ''%s'',' ...
                  ' ''%s'' or ''%s''.\n'], ...
                  'LM','SM','LA','SA','BE','LM','SM','LR','SR','LI','SI')
            end
            eigs_sigma = upper(eigs_sigma);
            if strcmp(eigs_sigma,'SM')
               % eigs('SM') => ARPACK which='LM', sigma=0
               whch = 'LM';
            else
               % eigs(string), where string~='SM' => ARPACK which=string, sigma=0
               whch = eigs_sigma;
            end
            sigma = zeros(classAB);
         else
            % eigs(scalar) => ARPACK which='LM', sigma=scalar
            if ~isfloat(eigs_sigma) || ~isscalar(eigs_sigma)
               error('MATLAB:eigs:EigenvalueShiftNonScalar',...
                  'Eigenvalue shift sigma must be a scalar.')
            end
            sigma = eigs_sigma;
            if issparse(sigma)
               sigma = full(sigma);
            end
            sigma = cast(sigma,classAB);
            isrealprob = isrealprob && isreal(sigma);
            whch = 'LM';
         end
      end

      % Process the input OPTS and derive some ARPACK values.
      % ARPACK's minimum tolerance is eps/2 ([S/D]LAMCH's EPS)
      tol = eps(classAB);
      maxit = [];
      % code by LVF:
      maxtime =[];
      % end code by LVF
      p = [];
      % Always use resid as the start vector, whether it is OPTS.v0 or
      % randomly generated within eigs.  We default resid to empty here.
      % If the user does not initialize it, we provide a random residual
      % below.
      info = int32(1);
      resid = []; 
      eigs_display = 1;
      cholB = false; % do we have B or its Cholesky factor?
      permB = []; % if cholB, is it chol(B), or chol(B(permB,permB))?
      if (nargin >= (6-argOffset))
         opts = varargin{6-argOffset};
         if ~isa(opts,'struct')
            error('MATLAB:eigs:OptionsNotStructure',...
               'Options argument must be a structure.')
         end
         if isfield(opts,'issym') && ~Amatrix
            issymA = opts.issym;
            if (issymA ~= false) && (issymA ~= true)
               error('MATLAB:eigs:InvalidOptsIssym', ...
                  'opts.issym must be true or false.')
            end
         end
         if isfield(opts,'isreal') && ~Amatrix
            if (opts.isreal ~= false) && (opts.isreal ~= true)
               error('MATLAB:eigs:InvalidOptsIsreal', ...
                  'opts.isreal must be true or false.')
            end
            isrealprob = isrealprob && opts.isreal;
         end
         if ~isempty(B) && (isfield(opts,'cholB') || isfield(opts,'permB'))
            if isfield(opts,'cholB')
               cholB = opts.cholB;
               if (cholB ~= false) && (cholB ~= true)
                  error('MATLAB:eigs:InvalidOptsCholB', ...
                     'opts.cholB must be true or false.')
               end
               if isfield(opts,'permB')
                  if issparse(B) && cholB
                     permB = opts.permB;
                     if ~isvector(permB) || ~isequal(sort(permB(:)),(1:n)')
                        error('MATLAB:eigs:InvalidOptsPermB',...
                           'opts.permB must be a permutation of 1:n.')
                     end
                  else
                     warning('MATLAB:eigs:IgnoredOptionPermB', ...
                        ['Ignoring opts.permB since B is not its sparse' ...
                        ' Cholesky factor.'])
                  end
               end
            end
         end
         if isfield(opts,'tol')
            if ~isfloat(tol) || ~isscalar(opts.tol) || ~isreal(opts.tol) || (opts.tol<=0)
               error('MATLAB:eigs:InvalidOptsTol',...
                  ['Convergence tolerance opts.tol must be a strictly' ...
                  ' positive real scalar.'])
            end
            tol = cast(full(opts.tol),classAB);
         end
         if isfield(opts,'p')
            p = opts.p;
            pstr = ['Number of basis vectors opts.p must be a positive' ...
               ' integer <= n.'];
            if ~isnumeric(p) || ~isscalar(p) || ~isreal(p) || (p<=0) || (p>n)
               error('MATLAB:eigs:InvalidOptsP', pstr)
            end
            if issparse(p)
               p = full(p);
            end
            if (round(p) ~= p)
               warning('MATLAB:eigs:NonIntegerVecQty',['%s\n         ' ...
                  'Rounding number of basis vectors.'],pstr)
               p = round(p);
            end
         end
         if isfield(opts,'maxit')
            maxit = opts.maxit;
            str = ['Maximum number of iterations opts.maxit must be' ...
               ' a positive integer.'];
            if ~isnumeric(maxit) || ~isscalar(maxit) || ~isreal(maxit) || (maxit<=0)
               error('MATLAB:eigs:OptsMaxitNotPosInt', str)
            end
            if issparse(maxit)
               maxit = full(maxit);
            end
            if (round(maxit) ~= maxit)
               warning('MATLAB:eigs:NonIntegerIterationQty',['%s\n         ' ...
                  'Rounding number of iterations.'],str)
               maxit = round(maxit);
            end
         end
         if isfield(opts,'v0')
            if ~isfloat(opts.v0) || ~isequal(size(opts.v0),[n,1])
               error('MATLAB:eigs:WrongSizeOptsV0',...
                  'Start vector opts.v0 must be n-by-1.')
            end
            if isrealprob
               if ~isreal(opts.v0)
                  error('MATLAB:eigs:NotRealOptsV0',...
                     'Start vector opts.v0 must be real for real problems.')
               end
               resid(1:n,1) = full(opts.v0);
            else
               resid(2:2:2*n,1) = full(imag(opts.v0));
               resid(1:2:(2*n-1),1) = full(real(opts.v0));
            end
         end
         if isfield(opts,'disp')
            eigs_display = opts.disp;
            dispstr = 'Diagnostic level opts.disp must be an integer.';
            if ~isnumeric(eigs_display) || ~isscalar(eigs_display) || ...
                  ~isreal(eigs_display) || (eigs_display<0)
               error('MATLAB:eigs:NonIntegerDiagnosticLevel', dispstr)
            end
            if (round(eigs_display) ~= eigs_display)
               warning('MATLAB:eigs:NonIntegerDiagnosticLevel', ...
                  '%s\n         Rounding diagnostic level.',dispstr)
               eigs_display = round(eigs_display);
            end
         end
         if isfield(opts,'cheb')
            error('MATLAB:eigs:ObsoleteOptionCheb', ...
               'Polynomial acceleration opts.cheb is an obsolete option.');
         end
         if isfield(opts,'stagtol')
            error('MATLAB:eigs:ObsoleteOptionStagtol', ...
               'Stagnation tolerance opts.stagtol is an obsolete option.');
         end
         % start code by LVF:
         if isfield(opts,'maxtime')
            maxtime = opts.maxtime;
            str = ['Maximum time opts.maxtime must be' ...
               ' a positive real number.'];
            if ~isnumeric(maxtime) || ~isscalar(maxtime) ||  ...
                    ~isreal(maxtime) || (maxtime<=0)
               error('MATLAB:eigs:OptsMaxTimeNotPos', str)
            end
            if issparse(maxtime)
               maxtime = full(maxtime);
            end            
         end
         % end code by LVF
      end
      if (isempty(resid))
        if isrealprob
           resid = cast(rand(n,1),classAB);
        else
           resid = cast(rand(2*n,1),classAB);
        end
      end

      afunNargs = zeros(1,0);
      if ~Amatrix
         % The trailing parameters for afun start at varargin{7-argOffset}
         % in eigs(afun,n,[B],k,sigma,opts,P1,P2,...). If there are no
         % trailing parameters in eigs, then afunNargs is a 1-by-0 empty
         % and no trailing parameters are passed to afun(x)
         afunNargs = 7-argOffset:nargin;
      end

      % Now that OPTS has been processed, do final error checking and
      % assign ARPACK variables

      % Extra check on input B
      if ~isempty(B)
         % B must be symmetric (Hermitian) positive (semi-)definite
         if cholB
            if ~isequal(triu(B),B)
               error('MATLAB:eigs:BsizeMismatchAorNotSPDorNotChol', Bstr)
            end
         else
            if ~ishermitian(B)
               error('MATLAB:eigs:BsizeMismatchAorNotSPDorNotChol', Bstr)
            end
         end
      end

      % Extra check on input K
      % We fall back on using the full EIG code if K is too large.
      useeig = false;
      if isrealprob && issymA
         knstr = sprintf(['For real symmetric problems, must have' ...
            ' number of eigenvalues k < n.\n']);
      else
         knstr = sprintf(['For nonsymmetric and complex problems,' ...
            ' must have number of eigenvalues k < n-1.\n']);
      end
      if isempty(B)
         knstr = [knstr 'Using eig(full(A)) instead.'];
      else
         knstr = [knstr 'Using eig(full(A),full(B)) instead.'];
      end
      if (k == 0)
         useeig = true;
      end
      if isrealprob && issymA
         if (k > n-1)
            if (n >= 6)
               warning('MATLAB:eigs:TooManyRequestedEigsForRealSym', ...
                  '%s',knstr)
            end
            useeig = true;
         end
      else
         if (k > n-2)
            if (n >= 7)
               warning('MATLAB:eigs:TooManyRequestedEigsForComplexNonsym', ...
                  '%s',knstr)
            end
            useeig = true;
         end
      end

      % Extra check on input SIGMA
      if isrealprob && issymA
         if ~isreal(sigma)
            error('MATLAB:eigs:ComplexShiftForRealSymProblem',...
               ['For real symmetric problems, eigenvalue shift sigma must' ...
               ' be real.'])
         end
      else
         if ~isrealprob && issymA && ~isreal(sigma)
            warning('MATLAB:eigs:ComplexShiftForHermitianProblem', ...
               ['Complex eigenvalue shift sigma on a Hermitian problem' ...
               ' (all real eigenvalues).'])
         end
      end
      if isrealprob && issymA
         if strcmp(whch,'LR')
            whch = 'LA';
            warning('MATLAB:eigs:SigmaChangedToLA', ...
               ['For real symmetric problems, sigma value ''LR''' ...
               ' (Largest Real) is now ''LA'' (Largest Algebraic).'])
         end
         if strcmp(whch,'SR')
            whch = 'SA';
            warning('MATLAB:eigs:SigmaChangedToSA', ...
               ['For real symmetric problems, sigma value ''SR''' ...
               ' (Smallest Real) is now ''SA'' (Smallest Algebraic).'])
         end
         if ~ismember(whch,{'LM', 'SM', 'LA', 'SA', 'BE'})
            error('MATLAB:eigs:EigenvalueRangeNotValid', ...
               [whchstr '\nFor real symmetric A, the' ...
               ' choices are ''%s'', ''%s'', ''%s'', ''%s'' or ''%s''.'], ...
               'LM','SM','LA','SA','BE');
         end
      else
         if strcmp(whch,'BE')
            warning('MATLAB:eigs:SigmaChangedToLM', ...
               ['Sigma value ''BE'' is now only available for real' ...
               ' symmetric problems.  Computing ''LM'' eigenvalues instead.'])
            whch = 'LM';
         end
         if ~ismember(whch,{'LM', 'SM', 'LR', 'SR', 'LI', 'SI'})
            error('MATLAB:eigs:EigenvalueRangeNotValid', ...
               [whchstr '\nFor non-symmetric or complex' ...
               ' A, the choices are ''%s'', ''%s'', ''%s'', ''%s'',' ...
               ' ''%s'' or ''%s''.\n'],'LM','SM','LR','SR','LI','SI');
         end
      end
      
      % The remainder of the error checking does not apply for the large
      % values of K that force us to use full EIG instead of ARPACK.
      if useeig
         return
      end

      % Extra check on input OPTS.p
      if isempty(p)
         if isrealprob && ~issymA
            p = min(max(2*k+1,20),n);
         else
            p = min(max(2*k,20),n);
         end
      else
         if isrealprob && issymA
            if (p <= k)
               error('MATLAB:eigs:InvalidOptsPforRealSymProb',...
                  ['For real symmetric problems, must have number of' ...
                  ' basis vectors opts.p > k.'])
            end
         else
            if (p <= k+1)
               error('MATLAB:eigs:InvalidOptsPforComplexOrNonSymProb',...
                  ['For nonsymmetric and complex problems, must have number of' ...
                  ' basis vectors opts.p > k+1.'])
            end
         end
      end

      % Extra check on input OPTS.maxit
      if isempty(maxit)
         maxit = max(300,ceil(2*n/max(p,1)));
      end

   end % checkInputs

%-------------------------------------------------------------------------%
   function fullEig(nOutputs)
      % Use EIG(FULL(A)) or EIG(FULL(A),FULL(B)) instead of ARPACK
      if ~isempty(B)
         B = Bmtimes(eye(n));
      end
      if isfloat(A)
         if issparse(A);
            A = full(A);
         end
      else
         % A is specified by a function.
         % Form the matrix A by applying the function.
         if ischar(eigs_sigma) && ~strcmp(eigs_sigma,'SM')
            % A is a function multiplying A*x
            AA = eye(n);
            for i = 1:n
               AA(:,i) = A(AA(:,i),varargin{afunNargs});
            end
            A = AA;
         else
            if (isfloat(eigs_sigma) && eigs_sigma == 0) || strcmp(eigs_sigma,'SM')
               % A is a function solving A\x
               invA = eye(n);
               for i = 1:n
                  invA(:,i) = A(invA(:,i),varargin{afunNargs});
               end
               A = eye(n) / invA;
            else
               % A is a function solving (A-sigma*B)\x
               % B may be [], indicating the identity matrix
               % U = (A-sigma*B)\sigma*B
               % => (A-sigma*B)*U = sigma*B
               % => A*U = sigma*B(U + eye(n))
               % => A = sigma*B(U + eye(n)) / U
               if isempty(B)
                  sB = eigs_sigma*eye(n);
               else
                  sB = eigs_sigma*B;
               end
               U = zeros(n,n);
               for i = 1:n
                  U(:,i) = A(sB(:,i),varargin{afunNargs});
               end
               A = sB*(U+eye(n)) / U;
            end
         end
      end

      if isempty(B)
         eigInputs = {A};
      else
         eigInputs = {A,B};
      end
      % Now with full floating point matrices A and B, use EIG:
      if (nOutputs <= 1)
         d = eig(eigInputs{:});
      else
         [V,D] = eig(eigInputs{:});
         d = diag(D);
      end

      % Grab the eigenvalues we want, based on sigma
      firstKindices = 1:k;
      lastKindices = n:-1:n-k+1;
      if ischar(eigs_sigma)
         switch eigs_sigma
            case 'LM'
               [ignore,ind] = sort(abs(d));
               range = lastKindices;
            case 'SM'
               [ignore,ind] = sort(abs(d));
               range = firstKindices;
            case 'LA'
               [ignore,ind] = sort(d);
               range = lastKindices;
            case 'SA'
               [ignore,ind] = sort(d);
               range = firstKindices;
            case 'LR'
               [ignore,ind] = sort(abs(real(d)));
               range = lastKindices;
            case 'SR'
               [ignore,ind] = sort(abs(real(d)));
               range = firstKindices;
            case 'LI'
               [ignore,ind] = sort(abs(imag(d)));
               range = lastKindices;
            case 'SI'
               [ignore,ind] = sort(abs(imag(d)));
               range = firstKindices;
            case 'BE'
               [ignore,ind] = sort(abs(d));
               range = [1:floor(k/2), n-ceil(k/2)+1:n];
            otherwise
               error('MATLAB:eigs:fullEigSigma','Unknown value of sigma');
         end
      else
         % sigma is a scalar
         [ignore,ind] = sort(abs(d-eigs_sigma));
         range = 1:k;
      end
      
      if (nOutputs <= 1)
         varargout{1} = d(ind(range));
      else
         varargout{1} = V(:,ind(range));
         varargout{2} = D(ind(range),ind(range));
         if (nOutputs == 3)
            % flag indicates "convergence"
            varargout{3} = 0;
         end
      end
      
   end % FULLEIG
   
%-------------------------------------------------------------------------%
   function [RB,RBT,perm] = CHOLfactorB
      % permB may be [] (from checkInputs) if the problem is not sparse
      % or if it was not passed in as opts.permB
      perm = permB;
      if cholB
         % CHOL(B) was passed in as B
         RB = B;
         RBT = B';
      else
         % CHOL(B) was not passed into EIGS
         if (mode == 1) && ~isempty(B)
            % Algorithm requires CHOL(B) to be computed
            if issparse(B)
               perm = symamd(B);
               [RB,pB] = chol(B(perm,perm));
            else
               [RB,pB] = chol(B);
            end
            if (pB == 0)
               RBT = RB';
            else
               error('MATLAB:eigs:BNotSPD', ...
                  'B is not symmetric positive definite.')
            end
         end
      end
   end % CHOLfactorB

%-------------------------------------------------------------------------%
   function [L,U,P,perm] = LUfactorAminusSigmaB
      % LU factor A-sigma*B, including a reordering perm if it is sparse
      if isempty(B)
         if issparse(A)
            AsB = A - sigma * speye(n);
         else
            AsB = A - sigma * eye(n);
         end
      else
         if cholB
            if issparse(B)
               AsB = A - sigma * Bmtimes(speye(n));
            else
               AsB = A - sigma * Bmtimes(eye(n));
            end
         else
            AsB = A - sigma * B;
         end
      end
      if issparse(AsB)
         [L,U,P,Q] = lu(AsB);
         [perm,ignore] = find(Q);
      else
         [L,U,P] = lu(AsB);
         perm = [];
      end
      % Warn if lu(A-sigma*B) is ill-conditioned
      % => sigma is close to an exact eigenvalue of (A,B)
      dU = diag(U);
      rcondestU = full(min(abs(dU)) / max(abs(dU)));
      if (rcondestU < eps)
         if isempty(B)
            ds = '(A-sigma*I)';
         else
            ds = '(A-sigma*B)';
         end
         warning('MATLAB:eigs:SigmaNearExactEig',...
            [ds ' has small reciprocal condition' ...
            ' estimate: %f\n' ...
            '         indicating that sigma is near an exact' ...
            ' eigenvalue.\n         The algorithm may not converge unless' ...
            ' you try a new value for sigma.\n'], ...
            rcondestU);
      end
   end % LUfactorAminusSigmaB

%-------------------------------------------------------------------------%
   function cols = checkIpntr
      % Check that ipntr returned from ARPACK refers to the start of a
      % column of workd.
      if ~isempty(B) && (mode == 3) && (ido == 1)
         inds = double(ipntr(1:3));
      else
         inds = double(ipntr(1:2));
      end
      [rows,cols] = ind2sub([n,3],inds);
      nonOneRows = find(rows~=1);
      if ~isempty(nonOneRows)
         error('MATLAB:eigs:ipntrMismatchWorkdColumn', ...
         ['One of ipntr(1:3) does not refer to the start' ...
            ' of a column of the %d-by-3 array workd.'],n)
      end
   end % checkIpntr

%-------------------------------------------------------------------------%
   function v = Amtimes(u)
      % Matrix-vector multiply v = A*u
      if Amatrix
         v = A * u;
      else % A is a function
         v = A(u,varargin{afunNargs});
         if isrealprob && ~isreal(v)
            error('MATLAB:eigs:complexFunction', ...
                  'AFUN is complex; set opts.isreal = false.');
         end
      end
   end

%-------------------------------------------------------------------------%
   function v = Bmtimes(u)
      % Matrix-vector multiply v = B*u
      if cholB % use B's cholesky factor and its transpose
         if ~isempty(permB)
            v(permB,:) = RBT * (RB * u(permB,:));
         else
            v = RBT * (RB * u);
         end
      else
         v = B * u;
      end
   end

%-------------------------------------------------------------------------%
   function v = RBsolve(u)
      % Solve v = RB\u for v
      if issparse(B)
         if ~isempty(permB)
            v(permB,:) = RB \ u;
         else
            v = RB \ u;
         end
      else
         RBopts.UT = true;
         v = linsolve(RB,u,RBopts);
      end
   end

%-------------------------------------------------------------------------%
   function v = RBTsolve(u)
      % Solve v = RB'\u for v
      if issparse(B)
         if ~isempty(permB)
            v = RBT \ u(permB,:);
         else
            v = RBT \ u;
         end
      else
         RBTopts.LT = true;
         v = linsolve(RBT,u,RBTopts);
      end
   end

%-------------------------------------------------------------------------%
   function v = AminusSigmaBsolve(u)
      % Solve v = (A-sigma*B)\u for v
      if Amatrix
         if ~isempty(permAsB)
            % use LU reordering permAsB
            v(permAsB,:) = U \ (L \ (P * u));
         else
            v = U \ (L \ (P * u));
         end
      else % A is a function
         v = A(u,varargin{afunNargs});
         if isrealprob && ~isreal(v)
            error('MATLAB:eigs:complexFunction', ...
                  'AFUN is complex; set opts.isreal = false.');
         end
      end
   end % AminusSigmaBsolve

%-------------------------------------------------------------------------%
   function displayRitzValues
      % Display a few Ritz values at the current iteration
      iter = double(ipntr(15));
      if (iter > eigs_iter) && (ido ~= 99)
         eigs_iter = iter;
         ds = sprintf(['Iteration %d: a few Ritz values of the' ...
            ' %d-by-%d matrix:'],iter,p,p);
         disp(ds)
         if isrealprob
            if issymA
               dispvec = workl(double(ipntr(6))+(0:p-1));
               if strcmp(whch,'BE')
                  % roughly k Large eigenvalues and k Small eigenvalues
                  disp(dispvec(max(end-2*k+1,1):end))
               else
                  % k eigenvalues
                  disp(dispvec(max(end-k+1,1):end))
               end
            else
               dispvec = complex(workl(double(ipntr(6))+(0:p-1)), ...
                  workl(double(ipntr(7))+(0:p-1)));
               % k+1 eigenvalues (keep complex conjugate pairs together)
               disp(dispvec(max(end-k,1):end))
            end
         else
            dispvec = complex(workl(2*double(ipntr(6))-1+(0:2:2*(p-1))), ...
               workl(2*double(ipntr(6))+(0:2:2*(p-1))));
            disp(dispvec(max(end-k+1,1):end))
         end
      end
   end

%-------------------------------------------------------------------------%
   function flag = processEUPDinfo(warnNonConvergence)
      % Process the info flag returned by the ARPACK routine **eupd
      flag = 0;
      if (info ~= 0)
         es = ['Error with ARPACK routine ' eupdfun ':\n'];
         switch double(info)
            case 2
               ss = sum(select);
               if (ss < k)
                  error('MATLAB:eigs:ARPACKroutineError02ssLTk', ...
                     [es 'The logical variable select was only set' ...
                     ' with %d 1''s instead of nconv=%d (k=%d).\n' ...
                     'Please report this to the ARPACK authors at' ...
                     ' arpack@caam.rice.edu.'], ...
                     ss,double(iparam(5)),k)
               else
                  error('MATLAB:eigs:ARPACKroutineError02', ...
                     [es 'The LAPACK reordering routine %strsen' ...
                     ' did not return all %d eigenvalues.'], ...
                     aupdfun(1),k);
               end
            case 1
               error('MATLAB:eigs:ARPACKroutineError01', ...
                  [es 'The Schur form could not be reordered by the' ...
                  ' LAPACK routine %strsen.\nPlease report this to the' ...
                  ' ARPACK authors at arpack@caam.rice.edu.'], ...
                  aupdfun(1))
            case -14
               error('MATLAB:eigs:ARPACKroutineErrorMinus14', ...
                  [es aupdfun ...
                  ' did not find any eigenvalues to sufficient accuracy.']);
            otherwise
               error('MATLAB:eigs:ARPACKroutineError', ...
                  [es 'info = %d. Please consult the ARPACK Users''' ...
                  ' Guide for more information.'],full(info));
         end
      else
         nconv = double(iparam(5));
         if (nconv == 0)
            if (warnNonConvergence)
               warning('MATLAB:eigs:NoEigsConverged', ...
                  'None of the %d requested eigenvalues converged.',k)
            else
               flag = 1;
            end
         elseif (nconv < k)
            if (warnNonConvergence)
               warning('MATLAB:eigs:NotAllEigsConverged', ...
                  'Only %d of the %d requested eigenvalues converged.', ...
                  nconv,k)
            else
               flag = 1;
            end
         end
      end
   end % processEUPDinfo

%-------------------------------------------------------------------------%
   function printTimings
      % Print the time taken for each major stage of the EIGS algorithm
      if (mode == 1)
         innerstr = sprintf(['Compute A*X:' ...
            '                               %f\n'],cputms(3));
      elseif (mode == 3)
         if isempty(B)
            innerstr = sprintf(['Solve (A-SIGMA*I)*X=Y for X:' ...
               '               %f\n'],cputms(3));
         else
            innerstr = sprintf(['Solve (A-SIGMA*B)*X=B*Y for X:' ...
               '             %f\n'],cputms(3));
         end
      end
      if ((mode == 3) && (Amatrix))
         if isempty(B)
            prepstr = sprintf(['Pre-processing, including lu(A-sigma*I):' ...
               '   %f\n'],cputms(1));
         else
            prepstr = sprintf(['Pre-processing, including lu(A-sigma*B):' ...
               '   %f\n'],cputms(1));
         end
      else
         prepstr = sprintf(['Pre-processing:' ...
            '                            %f\n'],cputms(1));
      end
      sstr = sprintf('***********CPU Timing Results in seconds***********');
      ds = sprintf(['\n' sstr '\n' ...
         prepstr ...
         'ARPACK''s %s:                           %f\n' ...
         innerstr ...
         'Post-processing with ARPACK''s %s:      %f\n' ...
         '***************************************************\n' ...
         'Total:                                     %f\n' ...
         sstr '\n'], ...
         aupdfun,cputms(2),eupdfun,cputms(4),cputms(5));
      disp(ds)
   end % printTimings

%-------------------------------------------------------------------------%
% End of nested functions
%-------------------------------------------------------------------------%

end % EIGS

%-------------------------------------------------------------------------%
% Subfunctions
%-------------------------------------------------------------------------%
function tf = ishermitian(A)
%ISHERMITIAN
tf = isequal(A,A');
end % ishermititan
%-------------------------------------------------------------------------%
% End of subfunctions
%-------------------------------------------------------------------------%







