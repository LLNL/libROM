function [inv_Q] = gnat(Q, m_used, nsr)
% Input:      Q: n by m matrix with orthonormal columns
%        m_used: number of basis_vectors to use
%           nsr: number of samples required

    phi = [];
    used = [];
    U = Q(:,1);
    Q_sampled = [];
    m = min(m_used, size(Q,2));
    if nsr > 0
        ns = nsr;
    else
        ns = m;
    end
    n = size(Q,1);
    ns_mod_nr = mod(ns,m);
    if 0 < ns_mod_nr
        nsi = idivide(int32(ns), int32(m), 'floor') + 1;
    else
        nsi = idivide(int32(ns), int32(m), 'floor');
    end
    P = [];
    for i = 1:nsi
        s_row = -1;
        s_row_val = Inf;
        for j = 1:n
            if ~any(used(:) == j)
                if s_row == -1 || s_row_val < abs(Q(j,1))
                    s_row = j;
                    s_row_val = abs(Q(j,1));
                end
            end
        end
        used(end + 1) = s_row;
        phi(end + 1) = s_row;
        newPcol = zeros(n,1);
        newPcol(phi(end)) = 1;
        P = [P newPcol];
        Q_sampled = [Q_sampled;Q(phi(end),:)];
    end
    for l = 2:m
        M = transpose(P)*U;
        inv_M = pinv(M);
        RHS = transpose(P) * Q(:, l);
        c = inv_M * RHS;
        r = Q(:,l) - U*c;
        U = [U Q(:,l)];
        if l - 1 < ns_mod_nr
            nsi = idivide(int32(ns), int32(m), 'floor') + 1;
        else
            nsi = idivide(int32(ns), int32(m), 'floor');
        end
        for i = 1:nsi
            s_row = -1;
            s_row_val = Inf;
            for j = 1:n
                if ~any(used(:) == j)
                    if s_row == -1 || s_row_val < abs(r(j,1))
                        s_row = j;
                        s_row_val = abs(r(j,1));
                    end
                end
            end
            used(end + 1) = s_row;
            phi(end + 1) = s_row;
            newPcol = zeros(n,1);
            newPcol(phi(end)) = 1;
            P = [P newPcol];
            Q_sampled = [Q_sampled;Q(phi(end),:)];
        end
    end
    [~, phi_sort_order] = sort(phi);
    Q_sampled = Q_sampled(phi_sort_order,:);
    inv_Q = transpose(pinv(Q_sampled));
end
