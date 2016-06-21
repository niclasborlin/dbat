function pts=ddegetallpts(ch,ids)
%pts=ddegetallpts(ch,ids)
%ch  - dde channel.
%ids - ids of points to try to get.
%pts - [id,x,y,z,p,c], where
%      p - processed (1), or approximated (0),
%      c - ctrl pt.

% Preallocate to avoid memory fragmentation.
pts=zeros(length(ids),6);
found=0;
for i=1:length(ids)
    cmd=sprintf('GetPoint %d',full(ids(i)));
    [ok,num,str]=ddecmd(ch,cmd);
    if (ok)
        found=found+1;
        pts(found,:)=[num{2},num{5},num{6},num{7},num{3},num{4}];
    end
end
% Free unused memory.
pts(found+1:end,:)=[];
