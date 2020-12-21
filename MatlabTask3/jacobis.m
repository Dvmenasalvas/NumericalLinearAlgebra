function [x,k,resvec,DD,ID,JD,D,Ab] = jacobis(AA,IA,JA,b,x,tol,kmax)

%This function is an implementation of Jacobi's
%iterative method for solving Ax=b.
%The input vector x is the initial guess, while tol
%is the tolerance used for stopping the iteration; kmax is the
%maximum number of iterations to be run. The output is the kth
%approximation to the exact solution, while resvec is a vector
%containing the 2-norms of the residual vectors.
n=length(b);
if ~exist('x'), x=zeros(n,1);end
if ~exist('tol'), tol=1e-6;end
if ~exist('kmax'), kmax=1e4;end

Ab=matvecs(AA,IA,JA,b);%ignore this line: used for testing
r=b-matvecs(AA,IA,JA,x);
res0=norm(r);res=norm(r);
resvec=res/res0;res=1;k=0;
%extract a vector D containing the main diagonal of matrix A
[DD,ID,JD]=diags(AA,IA,JA);
D=DD;
k = 1;
while k<length(DD)
   if ID(k)+1 ~= ID(k+1)
       zer = zeros(ID(k+1) - (ID(k)+1),1);
       D = [D(1:k); zer; D(k+1)];
       k = k + 1 + (ID(k+1) - (ID(k)+1));
   else, k = k + 1;
   end
end
while res>tol & k<kmax
    k=k+1;
    z=r./D;
    x=x+z;
    r=b-matvecs(AA,IA,JA,x);
    res=norm(r)/res0;
    resvec=[resvec res];
    fprintf('%2.9f\t  %3d\n',resvec(end),k)  
end

