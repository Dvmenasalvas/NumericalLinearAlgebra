function [DD,ID,JD]=diags(AA,IA,JA)

%Jacobi splitting:
%Input: sparse matrix A in coordinate format (AA, IA, JA)
%Output: sparse diagonal matrix D in coordinate format (DD, ID, JD)

DD=AA;
ID=IA;
JD=JA;
k = 1;
while k<=length(DD)
   if ID(k) ~= JD(k)
       DD(k) = [];
       ID(k) = [];
       JD(k) = [];
   else, k = k + 1;
   end
end
ID(end) = length(DD);
JD(end) = length(DD);