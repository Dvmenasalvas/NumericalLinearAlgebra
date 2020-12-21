%1: Exponentially-stretched grid
m=12;n=100;
beta=log(m);
%Edit the line below to generate the vector v
x = linspace(0,1,n)';
v = m*(1.-exp(-beta*x))/(m-1);

%----------------------------------------------------------------------
%2: Antidiagonal of a matrix
rng(1,'twister');n=200;
M=randi(10,n,n);
for i=1:n
    d(i,1)=M(i,n+1-i);
end
%----------------------------------------------------------------------
%3 Matrix-matrix product
m=12;n=15;q=9;
rng(1,'twister');A=randi(10,m,n);
rng(1,'twister');B=randi(10,n,q);
C=zeros(m,q);
for k=1:n
    %Edit the line below to generate the correct matrix Q
    Q=A(:,k)*B(k,:);
    C=C+Q;
end
%-----------------------------------------------------------------------
%4 Iterating with Gauss matrices
n=5;T=toeplitz(n:-1:1);
for k=1:n-1
    Lk=eye(n);
    Lk(k+1:n,k)= -T(k+1:n,k)/T(k,k);
    T=Lk*T;
end
%-----------------------------------------------------------------------
%5 LU for tridiagonal matrices
J=[];I=[];%ignore this line: used for checking your implementation
A=randi(9,10,10);A=triu(tril(A,1),-1);
n=length(A);
for k=1:n-1
    j=k+1;
    %There is no way to check wether A(k,k) is nonzero, i.e. if
    %LU descomposition is possible. As that is not mention in the question,
    %I decided not to check it, but that could be a improvement to
    %the programme.
    A(j,k)=A(j,k)/A(k,k);
        J=[J j];%ignore this line: used for checking your implementation
    i=k+1;
    A(j,i)=A(j,i)-A(j,k)*A(k,i);
            I=[I i];%ignore this line: used for checking your implementation
end
ind=length(I)*length(J);%ignore this line: used for checking your implementation
%--------DO NOT WRITE BELOW THIS LINE-----------------------------------------
