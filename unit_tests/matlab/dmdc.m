function dmdc(X, Y, r, t, dt, varargin)
% Input:      X: n by m snapshot matrix
%             Y: n by (m-1) control snapshot matrix
%             r: number of basis_vectors to use
%             t: The time to predict
%            dt: delta time
%             B: (optional) known control matrix

    X_in = X(:,1:end-1);
    X_out = X(:,2:end);
    b0 = X(:,1);

    if nargin == 6 % known control matrix
        B = varargin{1};
        X_out = X_out - B*Y;
    else
        X_in = [X_in; Y];
    end

    [U, S, V] = svd(X_in, 'econ');
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);

    m = size(X,1);
    U1 = U(1:m,:);
    if nargin == 6
        U_out = U1;
        Atilde = U_out'*X_out*V*inv(S);
        [W,ev] = eig(Atilde);
        Phi = X_out*V*inv(S)*W;
        Btilde = U_out'*B;
    else
        [U_out, ~, ~] = svd(X_out, 'econ');
        U_out = U_out(:,1:r);
        Atilde = U_out'*X_out*V*inv(S)*U1'*U_out;
        [W,ev] = eig(Atilde);
        Phi = X_out*V*inv(S)*U1'*U_out*W;
        U2 = U(m+1:end,:);
        Btilde = U_out'*X_out*V*inv(S)*U2';
    end

    n = size(X_out,2);
    pred = U_out*W*sum((diag(ev).^(n:-1:0)).*(W\[U_out'*b0, Btilde*Y]),2);

    norm(real(pred-X(:,end)))/norm(real(X(:,end)))
end
