function [inv_Q] = deim(Q)
% Input: Q: n by m matrix with orthonormal columns

    phi = [];
    [rho, phi(end + 1)] = max(abs(Q(:,1)));
    U = Q(:,1);
    m = size(Q,2);
    n = size(Q,1);
    P = zeros(n,1);
    P(phi(end)) = 1;
    Q_sampled = Q(phi(end),:);
    for l = 2:m
        M = transpose(P)*U;
        
        inv_M = inv(M);
        RHS = transpose(P) * Q(:, l);
        c = inv_M * RHS;
        r = Q(:,l) - U*c;
        [rho, phi(end + 1)] = max(abs(r));
        U = [U Q(:,l)];
        newPcol = zeros(n,1);
        newPcol(phi(end)) = 1;
        P = [P newPcol];
        Q_sampled = [Q_sampled;Q(phi(end),:)];
    end
    [phi, phi_sort_order] = sort(phi);
    Q_sampled = Q_sampled(phi_sort_order,:);
    inv_Q = inv(Q_sampled);
end
