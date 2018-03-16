function fail=compare_jacobian_structs(caller,A,N,absThres,relThres,verbose)
%COMPARE_JACOBIAN_STRUCTS Compare Jacobian structs.
%
%   FAIL=COMPARE_JACOBIAN_STRUCTS(CALLER,A,N,ABSTHRES,RELTHRES,VERBOSE)
%   Compare all fields of the structs A and N holding analytical and
%   numerical Jacobians, respectively. Return TRUE and output warnings
%   if any absolute or relative errors exceed given thresholds. If
%   VERBOSE, always output error levels.

fail=false;
fn=union(fieldnames(A),fieldnames(N));

for i=1:length(fn)
    if ~isfield(A,fn{i})
        fail=true;
        warning([caller,': Bug: %s not a field in analytical struct'],fn{i});
    elseif ~isfield(N,fn{i})
        fail=true;
        warning([caller,': Bug: %s not a field in numerical struct'],fn{i});
    elseif any(size(A.(fn{i}))~=size(N.(fn{i})))
        fail=true;
        warning([caller,': Size mismatch']);
    else
        absErr=abserr(A.(fn{i}),N.(fn{i}));
        relErr=relerr(A.(fn{i}),N.(fn{i}));
        if verbose
            fprintf('%s: %s absErr = %g, relErr= %g.\n',caller,fn{i},absErr,...
                    relErr);
        end
        if absErr>absThres
            fail=true;
            warning('%s: %s absolute error (%g) above threshold (%g)',...
                    caller,fn{i},absErr,absThres);
        elseif relErr>relThres
            fail=true;
            warning('%s: %s relative error (%g) above threshold (%g)',...
                    caller,fn{i},relErr,relThres);
        end
    end
end
