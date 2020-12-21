function [H,Q,c,s,B]=givensqr(H)
%Edit the lines below to provide your personal information
%Name: [Daniel Valverde Menasalvas]
%Student ID: [DXV925]
%This file computes the QR factorisation
%of a nxm Hessenberg matrix using Givens rotations.
%Please ignore output variables c,s,B: these are used for testing only.
[n,m]=size(H);
Q=eye(n);
for k=1:m
    a=H(k,k);b=H(k+1,k);
    [c,s]=givens(a,b);
    H=gmatmat(H,c,s,k);
    Gk = eye(n);
    Gk(k:k+1,k:k+1) = [c,-s;s,c];
    Q=Gk*Q;
end
Q = Q';
%Ignore the next line: used for testing
[c,s]=givens(1,1);B=gmatmat(fliplr(vander(1:3)),1/sqrt(2),-1/sqrt(2),1);
end
%---------------------------------

