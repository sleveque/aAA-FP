function [x_new, z_new, u_new, iter, res, norm_res] = ADMM(A, B, b, rho, z0, u0, tol, maxit, x_update, z_update, u_update)

    z_old = z0;
    u_old = u0;

    x_new = x_update(z_old, u_old);
    z_new = z_update(x_new, u_old);
    u_new = u_update(x_new, z_new, u_old);

    res = norm([z_new - z_old; u_new - u_old]);
    norm_res = zeros(1, maxit);
    norm_res(1) = res;
    fprintf('initial residual norm: %e\n', res)

    z_old = z_new;
    u_old = u_new;

    for iter = 2 : maxit
        if res < tol
            break
        end

        x_new = x_update(z_old, u_old);
        z_new = z_update(x_new, u_old);
        u_new = u_update(x_new, z_new, u_old);

        res = norm([z_new - z_old; u_new - u_old]);
        fprintf('iter % i, residual norm: %e\n', iter, res)

        norm_res(iter) = res;

        z_old = z_new;
        u_old = u_new;
    end

    norm_res = norm_res(1 : iter);

end
