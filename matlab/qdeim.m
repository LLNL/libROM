function [inv_U] = qdeim(U, m)
    [~, ~, p] = qr(U', 'vector');
    p = p(1:size(U, 2))';
    for i = length(p) + 1:m
        [~, S, W] = svd(U(p, :), 0);
        g = S(end - 1, end - 1).^2 - S(end, end)^2;
        Ub = W'*U';
        r = g + sum(Ub.^2,1);
        r = r - sqrt((g+sum(Ub.^2, 1)).^2 - 4 * g * Ub(end, :).^2);
        [~, I] = sort(r, 'descend');
        e = 1;
        while any(I(e) == p)
            e = e + 1;
        end
        p(end + 1) = I(e);
    end
    [p, ~] = sort(p);
    U_sampled = [];
    for i = 1:size(p,1)
        U_sampled = [U_sampled; U(p(i,1),:)];
    end
    if m == size(U,2)
        inv_U = pinv(U_sampled);
    else
        inv_U = transpose(pinv(U_sampled));
    end
end
