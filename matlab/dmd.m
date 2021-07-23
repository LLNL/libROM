function dmd(X, m_used)
% Input:      X: n by m snapshot matrix
%        m_used: number of basis_vectors to use

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);
    b0 = X(:,1);

    [U, S, V] = svd(X1, 'econ');
    U = U(:,1:m_used);
    U(1:5,1:5);
    V(1:5,1:5);
    S;
    inv(S)
    S = S(1:m_used,1:m_used);
    V = V(:,1:m_used);
    Atilde = U'*X2*V*inv(S);
    [W,eigs] = eig(Atilde);
    eigs;
    Phi = X2*V*inv(S)*W;
    v = Phi*eigs^100*pinv(Phi)*b0;
    norm(real(v-X(:,end)))/2340

end
