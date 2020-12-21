function B=gmatmat(A,c,s,k)

%This file computes the product B=G*A
%where G is the Givens rotation defined by the input c,s.
%The product is performed efficiently by modifying 
%only rows k and k+1 in the matrix A 
%via a multiplication by the 2x2 matrix G.

G=[c,-s;s,c];
B=A;
B(k:k+1,:)= G*B(k:k+1,:);
end
%---------------------------------