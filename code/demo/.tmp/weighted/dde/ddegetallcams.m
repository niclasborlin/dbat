function cams=ddegetallcams(ch,camIds)
%cams=ddegetallcams(ch,ids)
%ch   - dde channel.
%ids  - ids of photos to try to get.
%cams - [id,x,y,z,omega,phi,kappa] with angles in radians.

cams=zeros(0,7);
for i=1:length(camIds)
    cmd=sprintf('GetPhotoStation %d 0',camIds(i));
    [ok,num,str]=ddecmd(ch,cmd);
    if (ok)
        cams(end+1,:)=[camIds(i),num{2},num{3},num{4},num{5},num{6},num{7}];
    end
end

