function [S, M] = qdeim(U)
% Input:  U: n by m matrix with orthonormal columns
% Output: S: selection of m row indices with guaranteed upper bound
%            norm(inv(U(S,:))) <= sqrt(n - m + 1) * O(2 ^ m)
%         M: the matrix U * inv(U(S, :));
%            The Q-DEIM projection of an n-by-1 vector f is M*f(S)
    [n, m] = size(U);
    if nargout == 1
        [~, ~, P] = qr(U', 'vector'); S = P(1:m);
    else
        [Q, R, P] = qr(U', 'vector'); S = P(1:m);
        M = [eye(m); (R(:, 1:m) \ R(:, (m + 1):n))'];
        Pinverse(P) = 1:n; M = M(Pinverse,:);
    end
end