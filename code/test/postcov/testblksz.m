load tmp

UA=inv(LA);

LBUA=LB*UA;

% Try different block sizes to compute UB and sum(UB.^2,1)

blks=ceil(logspace(log10(50000),log10(size(UA,2)),5));
times=nan(size(blks));

for i=1:length(blks)
    blk=blks(i)

    start=now;
    dp=sum(UA.^2,1)';

    for base=1:blk:size(UA,2)
        % Columns in this block.
        ix=base:min(base+blk-1,size(UA,2));
        
        LBUAblk=full(LBUA(:,ix));
        UB=-LC\LBUAblk;
        dp(ix)=dp(ix)+sum(UB.^2,1)';
    end
    
    times(i)=(now-start)*86400
end

