function test_all

tests={'scale2','scale3','xlat2','xlat3','lin3','lin2','pinhole', ...
       'affine2mat','affine2','lens_rad2','power_vec','rad_scale', ...
       'tang_scale','brown_rad','brown_tang','brown_dist_abs', ...
       'brown_dist_rel','xformpt2cam'};

%for i=length(tests):-1:1
for i=1:length(tests)
    fprintf('Testing %s...',tests{i});
    fail=feval(tests{i},'selftest');
    if fail
        fprintf('FAIL.\n');
    else
        fprintf('OK.\n');
    end
end
