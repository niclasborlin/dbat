function test_all

tests={'xlat2','xlat3','scale2','scale3','lin3','lin2','pinhole', ...
       'affine2mat','affine2','lens_rad2','power_vec','rad_scale', ...
       'tang_scale','brown_rad','brown_tang','brown_dist_abs', ...
       'brown_dist_rel','eulerrotmat','world2cam','eulerpinhole','affscale2'};

for i=length(tests):-1:1
    fprintf('Testing %s...',tests{i});
    fail=feval(tests{i},'selftest');
    if fail
        fprintf('FAIL.\n');
    else
        fprintf('OK.\n');
    end
end
