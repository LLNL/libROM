function [Phi] = dmd(X, m_used)
% Input:      X: n by m snapshot matrix
%        m_used: number of basis_vectors to use

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);

    [U, S, V] = svd(Q, 'econ');
    U = U(:,1:m_used);
    S = S(1:m_used,1:m_used);
    V = V(:,1:m_used);
    Atilde = U'*X2*V*inv(S);
    [W,eigs] = eig(Atilde);
    Phi = X2*V*inv(S)*W;

end
