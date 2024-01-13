function [U, T] = jacoby(K, F, N, sigma, err, crit)
    D = spdiags(diag(K), 0, N - 1, N - 1) ^ -1;
    U = zeros(N - 1, 1);

    tic
    for i = 1:1:crit
        % Сдви
        delta = sigma * D * (K * U - F);
        if (max(abs(delta)) > err)
            U = U - delta;
        else
            break;
        end
    end

    T = toc;
end