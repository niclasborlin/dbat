function L=cholesky_element(A)

n=size(A,1);

L=zeros(n);

for j=1:n
    L(j,j)=sqrt(A(j,j)-sum(L(j,1:j-1).^2));
    for i=j+1:n
        L(i,j)=(A(i,j)-sum(L(i,1:j-1).*L(j,1:j-1)))/L(j,j);
    end
end
