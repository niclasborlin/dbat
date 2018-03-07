tests={'scale2','scale3','xlat2','xlat3','lin3','lin2','pinhole','affine2mat','affine2'};

for i=length(tests):-1:1
    fprintf('Testing %s...',tests{i});
    fail=feval(tests{i},'selftest');
    if fail
        fprintf('FAIL.\n');
    else
        fprintf('OK.\n');
    end
end
