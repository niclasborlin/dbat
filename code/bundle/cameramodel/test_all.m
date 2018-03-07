tests={'scale2','scale3','xlat2','xlat3'};

for i=1:length(tests)
    fprintf('Testing %s...',tests{i});
    fail=feval(tests{i},'selftest');
    if fail
        fprintf('FAIL.\n');
    else
        fprintf('OK.\n');
    end
end
