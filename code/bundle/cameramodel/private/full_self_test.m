function fail=full_self_test(caller,params,absThres,relThres,verbose,testParams)
%FULL_SELF_TEST Run full self-test on DBAT model functions.
%
%   FULL_SELF_TEST(CALLER,PARAMS,ABSTHRES,RELTHRES,VERBOSE) runs a
%   full self-test on the analytical and numerical Jacobians of
%   CALLER. The cell array PARAMS contain the basic parameters to
%   be used for the self-test. ABSTHRES and RELTHRES are thresholds
%   of when to signal errors. All combinations of extra logical
%   parameters to CALLER are tested. Returns TRUE if any test
%   fails, otherwise FALSE. If VERBOSE, output all tests.
%
%   FULL_SELF_TEST(CALLER,PARAMS,ABSTHRES,RELTHRES,VERBOSE,N) will
%   only run combinations of the N first logical parameters.
%
%SEE ALSO: COMPARE_JACOBIAN_STRUCTS.

if nargin<6, testParams=length(params); end

% Generate all combination of zero or more extra logical
% parameters.
computeJacobians={{}};
for l=1:testParams
    words=0:pow2(l)-1;
    for i=1:length(words)
        computeJacobians{end+1}=num2cell(logical(bitget(words(i),1:l))); %#ok<AGROW>
    end
end

fail=false;

for i=1:length(computeJacobians)
    compute=computeJacobians{i};
    v0=feval(caller,params{:});
    [v1,A,N]=feval(caller,params{:},compute{:});

    absErr=abserr(v0,v1);
    relErr=relerr(v0,v1);
    if verbose
        fprintf('%s: res absErr = %g, relErr= %g.\n',caller,absErr,relErr);
    end
    if absErr>absThres
        fail=true;
        warning('%s: res absolute error (%g) above threshold (%g)',...
                caller,absErr,absThres);
    elseif relErr>relThres
        fail=true;
        warning('%s: res relative error (%g) above threshold (%g)',...
                caller,relErr,relThres);
    end

    if fail, return, end
    
    fail=fail | compare_jacobian_structs(caller,A,N,absThres,relThres,verbose);

    if fail, return, end
end
