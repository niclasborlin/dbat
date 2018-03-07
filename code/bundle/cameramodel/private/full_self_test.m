function fail=full_self_test(caller,params,absThres,relThres,verbose)
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
%SEE ALSO: COMPARE_JACOBIAN_STRUCTS.

% Generate all combination of zero or more extra logical
% parameters.
computeJacobians={{}};
for l=1:length(params)
    words=0:pow2(l)-1;
    for i=1:length(words)
        computeJacobians{end+1}=num2cell(logical(bitget(words(i),1:l)));
    end
end

for i=1:length(computeJacobians)
    compute=computeJacobians{i};
    [v,A,N]=feval(caller,params{:},compute{:});
    fail=compare_jacobian_structs(caller,A,N,absThres,relThres,verbose);
    if fail, return, end
end
