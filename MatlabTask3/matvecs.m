function w=matvecs(AA,IA,JA,v)

%Sparse matrix-vector product c=A*v: the matrix is provided
%in coordinate format AA, IA, JA, with m=IA(end), n=JA(end).
%The input vector v is assumed to be full.
%The output vector w is assumed to be full.

ma=IA(end);na=JA(end);nb=length(v);
if na~=nb,error('Wrong sizes: cannot perform matrix-vector product');end

w=zeros(ma,1);
for k=1:length(AA)
   w(IA(k))=w(IA(k))+AA(k)*v(JA(k));%edit this line to provide the correct expression for w
end