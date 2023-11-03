function dmdc(X, Y, ef, t, dt, varargin)
% Input:      X: n by m snapshot matrix
%             Y: n by (m-1) control snapshot matrix
%            ef: energy fraction
%             t: The time to predict
%            dt: delta time
%             B: (optional) known control matrix

    X_in = X(:,1:end-1);
    X_out = X(:,2:end);
    b0 = X(:,1);

    if nargin == 5
        X_in = [X_in; Y];
    elseif nargin == 6
        B = varargin{1};
        X_out = X_out - B*Y;
    else
        error('Invalid input');
    end

    [U, S, V] = svd(X_in, 'econ');
    r = sum(cumsum(diag(S))/sum(diag(S)) < ef)
    
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);

    m = size(X,1);
    U1 = U(1:m,:);
    if nargin == 5
        U2 = U(m+1:end,:);
        [U_out, S_out, ~] = svd(X_out, 'econ');
        r_out = sum(cumsum(diag(S_out))/sum(diag(S_out)) < ef)
        U_out = U_out(:,1:r_out);
        Atilde = U_out'*X_out*V*inv(S)*U1'*U_out;
        Btilde = U_out'*X_out*V*inv(S)*U2';
        [W,ev] = eig(Atilde);
        %Phi = X_out*V*inv(S)*U1'*U_out*W;
    else
        U_out = U1;
        Atilde = U_out'*X_out*V*inv(S);
        Btilde = U_out'*B;
        [W,ev] = eig(Atilde);
        %Phi = X_out*V*inv(S)*W;
    end

    n = size(X_out,2);
    pred = U_out*W*sum((diag(ev).^(n:-1:0)).*(W\[U_out'*b0, Btilde*Y]),2);

    norm(real(pred-X(:,end)))/norm(real(X(:,end)))
end
