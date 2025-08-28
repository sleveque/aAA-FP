global dim_y dim_lambda x_sol y_sol lambda_sol

seed = 100000;
rng(seed)

% window for AA
w_aAR = 10;
% number of AA iterations
aa_iter = 10;
% number of FP iterations
fp_iter = 10;

tol = 1e-12;
maxit = 1000;

num_test = 1;

% TV
if num_test == 1
    % dimension of the problem
    n = 1000;
    mu = 10;

    xhat = rand(n,1);

    % building matrices and rhs
    A = [-eye(n - 1), zeros(n - 1, 1)] + [zeros(n - 1, 1), eye(n-1)];
    B = - eye(n - 1);
    b = zeros(n - 1, 1);

    % regularization parameter
    beta = 0.001 * norm(b, inf);

    % initial guess
    y0 = zeros(n - 1, 1);
    lambda0 = y0;
    dim_y = size(y0, 1);
    dim_lambda = size(lambda0, 1);

    % updates of solutions
    x_update = @(y, lambda)x_update_TV(y, lambda, A, mu, xhat);
    y_update = @(x, lambda)y_update_TV(x, lambda, A, mu, beta);
    lambda_update = @(x, y, lambda)lambda_update_TV(x, y, lambda, A, B, b);
end

% Lasso
if num_test == 2
    density = 0.01;

    % dimension of the problem
    m = 150;
    n = 300;

    lambda = 1;
    mu = 10;

    A_test = sprand(m, n, density);

    % building matrices and rhs
    A = eye(n);
    B = - eye(n);
    b = zeros(n, 1);

    xhat = rand(m, 1);

    % initial guess
    y0 = zeros(n, 1);
    lambda0 = zeros(n, 1);
    dim_y = size(y0, 1);
    dim_lambda = size(lambda0, 1);

    % updates of solutions
    x_update = @(z, u)x_update_Lasso(z, u, A_test, mu, xhat);
    y_update = @(x, u)y_update_Lasso(x, u, mu, lambda);
    lambda_update = @(x, z, u)lambda_update_Lasso(x, z, u);
end

% NNLS
if num_test == 3
    density = 0.01;

    % dimension of the problem
    m = 150;
    n = 300;

    mu = 2;

    A_test = sprand(m, n, density);

    % building matrices and rhs
    A = eye(n);
    B = - eye(n);
    b = zeros(n, 1);

    xhat = rand(m, 1);

    % initial guess
    y0 = zeros(n, 1);
    lambda0 = zeros(n, 1);
    dim_y = size(y0, 1);
    dim_lambda = size(lambda0, 1);

    % updates of solutions
    x_update = @(z, u)x_update_NNLS(z, u, A_test, mu, xhat);
    y_update = @(x, u)y_update_NNLS(x, u, mu);
    lambda_update = @(x, z, u)lambda_update_NNLS(x, z, u);
end



fprintf('\n')
% ADMM
fprintf('starting ADMM\n')
tic
[x_new, y_new, lambda_new, iter, res, norm_res] = ADMM(y0, lambda0, tol, maxit, x_update, y_update, lambda_update);
T_ADMM = toc;

fprintf('Total ADMM iter: %i\n', iter)
fprintf('ADMM res: %i\n', res)
fprintf('CPU time: %f\n', T_ADMM)


fprintf('\n')
% aAR-FP
fprintf('Starting aAR-FP.\n')

x0 = [y0; lambda0];

ADMMfun = @(z, u)ADMM(z, u, tol, 1, x_update, y_update, lambda_update);

fpiterfun = @(x)fixed_point(x, ADMMfun);

tic
[x, iter_aAR_FP, relres_aAR, norm_story_aAR_FP] = aAR_FP(x0, fp_iter, aa_iter, w_aAR, tol, maxit, fpiterfun);
T_aAR_FP=toc;

fprintf('Total aAR-FP iter: %i\n', iter_aAR_FP)
fprintf('AA res: %i\n', relres_aAR)
fprintf('CPU time: %f\n', T_aAR_FP)

if aa_iter == 0
    save("ADMM_TV.mat","norm_res")
elseif fp_iter == 0
    save("AA(10)_TV.mat","norm_story_aAR_FP")
    else
        save("aAA(10)[10]-FP[10]_TV.mat","norm_story_aAR_FP")
end


% TV
function x_new = x_update_TV(y_old, lambda_old, A, mu, xhat)

    sqrt_mu = sqrt(mu);

    M = [A; 1/sqrt_mu * eye(size(xhat,1))];
    rhs = [y_old - lambda_old; 1/sqrt_mu * xhat];

    x_new = M \ rhs;

end

function y_new = y_update_TV(x_new, lambda_old, A, mu, beta)

    h = A * x_new + lambda_old;

    y_new = sign(h) .* max(abs(h)-beta/mu, zeros(size(h)));
end

function lambda_new = lambda_update_TV(x_new, y_new, lambda_old, A, B, b)
    lambda_new = lambda_old + A * x_new + B * y_new - b;
end

% Lasso
function x_new = x_update_Lasso(y_old, lambda_old, A, mu, xhat)

    M = A' * A;
    M = M + mu * eye(size(M));
    rhs = A' * xhat + mu * (y_old - lambda_old);

    x_new = M \ rhs;

end

function y_new = y_update_Lasso(x_new, lambda_old, mu, beta)

    h = x_new + lambda_old;

    y_new = sign(h) .* max(abs(h)-beta/mu, zeros(size(h)));
end

function lambda_new = lambda_update_Lasso(x_new, y_new, lambda_old)
    lambda_new = lambda_old + x_new - y_new;
end


% NNLS
function x_new = x_update_NNLS(y_old, lambda_old, A, mu, xhat)

    sqrt_mu = sqrt(2 * mu);

    M = [A; 1/sqrt_mu * eye(size(A,2))];
    rhs = [1/sqrt_mu * xhat; y_old - lambda_old];

    x_new = M \ rhs;

end

function y_new = y_update_NNLS(x_new, lambda_old, mu)

    h = (1/mu) * (x_new + lambda_old);

    y_new = max(h, zeros(size(h)));
end

function lambda_new = lambda_update_NNLS(x_new, y_new, lambda_old)
    lambda_new = lambda_old + x_new - y_new;
end



function fp_x_new = fixed_point(fp_x_old, ADMMfun)
    global dim_y dim_lambda x_sol y_sol lambda_sol

    yold = fp_x_old(1 : dim_y);
    lambdaold = fp_x_old(dim_y + 1 : dim_y + dim_lambda);
    
    [x_new, y_new, lambda_new, ~, ~, ~] = ADMMfun(yold, lambdaold);
    x_sol = x_new;
    y_sol = y_new;
    lambda_sol = lambda_new;

    fp_x_new = [y_new; lambda_new];
end