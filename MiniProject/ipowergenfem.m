function [lam,v,res,k]=ipowergenfem(A,M,tol,kmax,v)
%Made by Daniel Valverde Menasalvas DXV925
%
%This function computes an approximation to the eigenvector 
%corresponding to the second smaller eigenvalue of A*v = lam*M*v,
%supposing the smallest eigenvalue of A is 0 with eigenvector ones(n,1).
%v is the starting vector and kmax is the max number of iterations.

n=length(A);
%As the eingenvectors of the GEVP are real, we can take v to be real
if ~exist('v'),v=rand(n+1,1);end
if ~exist('kmax'),kmax=n;end
if ~exist('tol'),tol=1e-6;end


alpha = alphaestimation(A,M);
vi = ones(n,1);
c = alpha*M*vi;
P = preconditioner(A,M);
A = [A , c ; c' , 0];
M = [M , zeros(n,1) ; zeros(1,n), 1];

v=v/norm(v);
k=0;res=1; 
while k<kmax & res>tol
    k=k+1;
    [v,flag] = minres(A,M*v,1e-6,kmax,P);
    v=v/norm(v);    
    lam=(v'*A*v)/(v'*M*v);
    res=norm(A*v-lam*M*v);
end
end


function P = preconditioner(A,M)
%Made by Daniel Valverde Menasalvas DXV925
%
% Generates a preconditioner for the linear system of our problem

P=A+M;
n = length(A);
P = [P , zeros(n,1) ; zeros(1,n), 1];
end

function alpha = alphaestimation(A,M)
%Made by Daniel Valverde Menasalvas DXV925
%
%This function computes a number greater or equal than the second smallest 
%eigenvalue of the GEVP: A*v = lam*M*v
%Note: We assume that all of the eigenvalues of the GEVP are possitive
%and that the eigenvector of the smallest eigenvalue is ones(n,1).

n = length(A);
x = [1 - n; ones(n-1,1)];
alpha = (x'*A*x)/(x'*M*x);
end