function [pp,f,K,P,a,u]=UnpackIO(IO,nK,nP)
%Unpack inner orientation vector.


pp=IO(1:2);
f=IO(3);
K=IO(3+(1:nK));
P=IO(3+nK+(1:nP));
a=IO(3+nK+nP+(1:2));
u=IO(3+nK+nP+6+(1:2));
