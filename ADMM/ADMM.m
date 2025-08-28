function [x_new, y_new, lambda_new, iter, res, norm_res] = ADMM(y0, lambda0, tol, maxit, x_update, y_update, lambda_update)

    y_old = y0;
    lambda_old = lambda0;

    x_new = x_update(y_old, lambda_old);
    y_new = y_update(x_new, lambda_old);
    lambda_new = lambda_update(x_new, y_new, lambda_old);

    res = norm([y_new - y_old; lambda_new - lambda_old]);
    norm_res = zeros(1, maxit);
    norm_res(1) = res;
    fprintf('initial residual norm: %e\n', res)

    y_old = y_new;
    lambda_old = lambda_new;

    for iter = 2 : maxit
        if res < tol
            break
        end

        x_new = x_update(y_old, lambda_old);
        y_new = y_update(x_new, lambda_old);
        lambda_new = lambda_update(x_new, y_new, lambda_old);

        res = norm([y_new - y_old; lambda_new - lambda_old]);
        fprintf('iter % i, residual norm: %e\n', iter, res)

        norm_res(iter) = res;

        y_old = y_new;
        lambda_old = lambda_new;
    end

    norm_res = norm_res(1 : iter);

end