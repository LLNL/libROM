function dmd(X, m_used, t, dt)
% Input:      X: n by m snapshot matrix
%        m_used: number of basis_vectors to use
%             t: The time to predict
%            dt: delta time

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);
    b0 = X(:,1);

    [U, S, V] = svd(X1, 'econ');
    U = U(:,1:m_used);
    S = S(1:m_used,1:m_used);
    V = V(:,1:m_used);
    Atilde = U'*X2*V*inv(S);
    [W,eigs] = eig(Atilde);
    Phi = X2*V*inv(S)*W;
    pred = Phi*eigs^(t/dt)*pinv(Phi)*b0;
    norm(real(pred-X(:,end)))/norm(real(X(:,end)))

end
