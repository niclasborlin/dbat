pause off
loadplotdemo
romabundledemo
romabundledemo_selfcal
romabundledemo_imagevariant
camcaldemo
camcaldemo2
camcaldemo_allmodels
camcaldemo_missing_obs
camcaldemo_1ray
camcaldemo_no_datum
ll={'c1','c2','s1','s2','s3','s4'};
for i=1:length(ll)
    for t=[false,true]
        ll{i}
        t
        prague2016_pm(ll{i},t)
    end
end
prague2016_ps('s5')
ps_postproc('')
stpierrebundledemo_ps
sxb_prior_eo(false)
sxb_prior_eo(true)
disp('Done.');
