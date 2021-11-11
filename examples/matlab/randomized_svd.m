function [U,Sigma,V] = randomized_svd(A,k)

%   "Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate
%    Matrix Decompositions" by N. Halko, P. G. Martinsson, and J. A. Tropp

[m,n] = size(A);
if nargin == 1
    k = (m>n)*n + (m<=n)*m;
end

if m > n
    Omega = randn(n,k);
    Y = A * Omega;
    [Q,~,~] = qr(Y,'vector');
    B = Q'*A;
    [Uhat,Sigma,V] = svd(B);
    U = Q*Uhat;
else
    Omega = randn(m,k);
    Y = A' * Omega;
    [Q,~,~] = qr(Y,'vector');
    B = Q'*A';
    [Uhat,Sigma,Vd] = svd(B);
    Ud = Q*Uhat;
    U = Vd;
    V = Ud;
end

end
