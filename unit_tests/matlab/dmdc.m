function dmdc(X, Y, r, t, dt, varargin)
% Input:      X: n by m snapshot matrix
%             Y: n by (m-1) control snapshot matrix
%             r: number of basis_vectors to use
%             t: The time to predict
%            dt: delta time
%             B: (optional) known control matrix

    X1 = X(:,1:end-1);
    X2 = X(:,2:end);
    b0 = X(:,1);

    if nargin == 6 % known control matrix
        B = varargin{1};
        X2 = X2 - B*Y;
    else
        X1 = [X1; Y];
    end

    [U, S, V] = svd(X1, 'econ');
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);

    if nargin == 6
        Atilde = U'*X2*V*inv(S);
        %[W,eigs] = eig(Atilde);
        %Phi = X2*V*inv(S)*W;
        Btilde = U'*B;
    else
        m = size(X,1);
        U1 = U(m,:);
        U2 = U(m+1:end,:);
        Atilde = U'*X2*V*inv(S)*U1'*U;
        %[W,eigs] = eig(Atilde);
        %Phi = X2*V*inv(S)*U1'*U*W;
        Btilde = U'*X2*V*inv(S)*U2';
    end

    pred = U1'*b0;
    n = size(X1);
    for k = 1:n
        pred = Atilde * pred + Btilde * Y(:,k);
    end
    pred = U1*pred;
    norm(real(pred-X(:,end)))/norm(real(X(:,end)))

end
