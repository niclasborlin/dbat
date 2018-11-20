function [pp,f,K,P,b]=unpackio(IO,nK,nP)
%Unpack inner orientation vector.

f=IO(1,:);
pp=IO(2:3,:);
b=IO(4:5,:);
K=IO(5+(1:nK));
P=IO(5+nK+(1:nP));
