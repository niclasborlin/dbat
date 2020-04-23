compile
load Ltmp
diary /tmp/x.txt
full(LA)
full(LB)
LC
C=icpc_mex(LA,LB,LC,3);
blk=2;
ix=(blk-1)*3+(1:3);
UA=inv(LA(ix,ix));
LBUA=LB(:,ix)*UA;
UB=LC\LBUA;
d11=sum(UA.^2,1)+sum(UB.^2,1);
dp=@(A,i,j)A(:,i:3:end).*A(:,j:3:end);
d12=sum(dp(UA,1,2),1)+sum(dp(UB,1,2),1);
d13=sum(dp(UA,1,3),1)+sum(dp(UB,1,3),1);
d23=sum(dp(UA,2,3),1)+sum(dp(UB,2,3),1);
diary off
