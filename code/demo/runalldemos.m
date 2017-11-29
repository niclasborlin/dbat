loadplotdemo
romabundledemo
camcaldemo
ll={'c1','c2','s1','s2','s3','s4'};
for i=1:length(ll)
    for t=[false,true]
        prague2016_pm(ll{i},t)
    end
end
prague2016_ps('s5')
ps_postproc('')
stpierrebundledemo_ps
disp('Done.');

      
