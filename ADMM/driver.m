global dim_z dim_u x_sol z_sol u_sol

seed = 100000;
rng(seed)

w_aAA = 0;
aa_iter = 0;
fp_iter = 1;

tol = 1e-12;
maxit = 1000;

% num_test = 1:   TV
% num_test = 2:   lasso
% num_test = 3:   NNLS
num_test = 1;

% TV
if num_test == 1
    n = 1000;
    rho = 10;
    y = rand(n,1);

    A = [-eye(n - 1), zeros(n - 1, 1)] + [zeros(n - 1, 1), eye(n-1)];
    B = - eye(n - 1);
    b = zeros(n - 1, 1);

    alpha = 0.001 * norm(b, inf);

    z0 = zeros(n - 1, 1);
    u0 = z0;
    dim_z = size(z0, 1);
    dim_u = size(u0, 1);

    x_update = @(z, u)x_update_TV(z, u, A, rho, y);
    z_update = @(x, u)z_update_TV(x, u, A, rho, alpha);
    u_update = @(x, z, u)u_update_TV(x, z, u, A, B, b);
end

% lasso
if num_test == 2
    density = 0.01;

    m = 150;
    n = 300;

    lambda = 1;
    rho = 10;

    A_test = sprand(m, n, density);

    A = eye(n);
    B = - eye(n);
    b = zeros(n, 1);

    y = rand(m, 1);

    z0 = zeros(n, 1);
    u0 = zeros(n, 1);
    dim_z = size(z0, 1);
    dim_u = size(u0, 1);

    x_update = @(z, u)x_update_lasso(z, u, A_test, rho, y);
    z_update = @(x, u)z_update_lasso(x, u, rho, lambda);
    u_update = @(x, z, u)u_update_lasso(x, z, u);
end

% NNLS
if num_test == 3
    density = 0.01;

    m = 150;
    n = 300;

    rho = 2;

    A_test = sprand(m, n, density);

    A = eye(n);
    B = - eye(n);
    b = zeros(n, 1);

    y = rand(m, 1);

    z0 = zeros(n, 1);
    u0 = zeros(n, 1);
    dim_z = size(z0, 1);
    dim_u = size(u0, 1);

    x_update = @(z, u)x_update_NNLS(z, u, A_test, rho, y);
    z_update = @(x, u)z_update_NNLS(x, u, rho);
    u_update = @(x, z, u)u_update_NNLS(x, z, u);
end



fprintf('\n')
% ADMM
fprintf('starting ADMM\n')
tic
[x_new, z_new, u_new, iter, res, norm_res] = ADMM(A, B, b, rho, z0, u0, tol, maxit, x_update, z_update, u_update);
T_ADMM = toc;

fprintf('Total ADMM iter: %i\n', iter)
fprintf('ADMM res: %i\n', res)
fprintf('CPU time: %f\n', T_ADMM)


fprintf('\n')
% aAA_FP
fprintf('Starting aAA_FP.\n')

x0 = [z0; u0];

ADMMfun = @(z, u)ADMM(A, B, b, rho, z, u, tol, 1, x_update, z_update, u_update);

fpiterfun = @(x)fixed_point(x, ADMMfun);

tic
[x, iter_aAA_FP, relres_aAA, norm_story_aAA_FP] = aAA_FP(x0, fp_iter, aa_iter, w_aAA, tol, maxit, fpiterfun);
T_aAA_FP=toc;

fprintf('Total aAA-FP iter: %i\n', iter_aAA_FP)
fprintf('AA res: %i\n', relres_aAA)
fprintf('CPU time: %f\n', T_aAA_FP)


% TV
function x_new = x_update_TV(z_old, u_old, A, rho, y)

    sqrt_rho = sqrt(rho);

    M = [A; 1/sqrt_rho * eye(size(y,1))];
    rhs = [z_old - u_old; 1/sqrt_rho * y];

    x_new = M \ rhs;

end

function z_new = z_update_TV(x_new, u_old, A, rho, alpha)

    h = A * x_new + u_old;

    z_new = sign(h) .* max(abs(h)-alpha/rho, zeros(size(h)));
end

function u_new = u_update_TV(x_new, z_new, u_old, A, B, b)
    u_new = u_old + A * x_new + B * z_new - b;
end

% lasso
function x_new = x_update_lasso(z_old, u_old, A, rho, y)

    M = A' * A;
    M = M + rho * eye(size(M));
    rhs = A' * y + rho * (z_old - u_old);

    x_new = M \ rhs;

end

function z_new = z_update_lasso(x_new, u_old, rho, lambda)

    h = x_new + u_old;

    z_new = sign(h) .* max(abs(h)-lambda/rho, zeros(size(h)));
end

function u_new = u_update_lasso(x_new, z_new, u_old)
    u_new = u_old + x_new - z_new;
end


% NNLS
function x_new = x_update_NNLS(z_old, u_old, A, rho, y)

    sqrt_rho = sqrt(2 * rho);

    M = [A; 1/sqrt_rho * eye(size(A,2))];
    rhs = [1/sqrt_rho * y; z_old - u_old];

    x_new = M \ rhs;

end

function z_new = z_update_NNLS(x_new, u_old, rho)

    h = (1/rho) * (x_new + u_old);

    z_new = max(h, zeros(size(h)));
end

function u_new = u_update_NNLS(x_new, z_new, u_old)
    u_new = u_old + x_new - z_new;
end


% fixed point iteration
function fp_x_new = fixed_point(fp_x_old, ADMMfun)
    global dim_z dim_u x_sol z_sol u_sol

    zold = fp_x_old(1 : dim_z);
    uold = fp_x_old(dim_z + 1 : dim_z + dim_u);
    
    [x_new, z_new, u_new, ~, ~, ~] = ADMMfun(zold, uold);
    x_sol = x_new;
    z_sol = z_new;
    u_sol = u_new;

    fp_x_new = [z_new; u_new];
end
