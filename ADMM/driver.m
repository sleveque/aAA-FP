global dim_y dim_lambda x_sol y_sol lambda_sol

seed = 100000;
rng(seed)

% window for AA
w_aAR = 3;
% number of AA iterations
aa_iter = 3;
% number of FP iterations
fp_iter = 3;

tol = 1e-12;
maxit = 1000;
num_initial_guess = 50;
norm_res_story_tot = zeros(maxit, num_initial_guess);
num_iter_story_tot = zeros(num_initial_guess, 1);
CPU_story_tot = zeros(num_initial_guess, 1);

num_test = 3;

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
    y0 = rand(n - 1, num_initial_guess);
%     y0 = zeros(n - 1, 1);
%     lambda0 = y0;
    lambda0 = rand(n - 1, num_initial_guess);
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
    y0 = rand(n, num_initial_guess);
%     y0 = zeros(n, 1);
%     lambda0 = zeros(n, 1);
    lambda0 = rand(n, num_initial_guess);
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
    y0 = rand(n, num_initial_guess);
%     y0 = zeros(n, 1);
    lambda0 = rand(n, num_initial_guess);
%     lambda0 = zeros(n, 1);
    dim_y = size(y0, 1);
    dim_lambda = size(lambda0, 1);

    % updates of solutions
    x_update = @(z, u)x_update_NNLS(z, u, A_test, mu, xhat);
    y_update = @(x, u)y_update_NNLS(x, u, mu);
    lambda_update = @(x, z, u)lambda_update_NNLS(x, z, u);
end


window_values = zeros(5,1);
AAit_values = window_values;
FPit_values = AAit_values;

window_values(1) = 10;
window_values(2) = 1;
window_values(3) = 3;
window_values(4) = 8;
window_values(5) = maxit;

AAit_values(1) = 15;
AAit_values(2) = 15;
AAit_values(3) = 10;
AAit_values(4) = 15;
AAit_values(5) = 10;

FPit_values(1) = 5;
FPit_values(2) = 5;
FPit_values(3) = 3;
FPit_values(4) = 3;
FPit_values(5) = 10;


fprintf('\n')
% ADMM
fprintf('starting ADMM\n')
for i = 1 : num_initial_guess
    y0i = y0(:, i);
    lambda0i = lambda0(:, i);
    tic
    [x_new, y_new, lambda_new, iter, res, norm_res] = ADMM(y0i, lambda0i, tol, maxit, x_update, y_update, lambda_update);
    T_ADMM = toc;
    
    norm_res_story_tot(1:iter, i) = norm_res;
    num_iter_story_tot(i) = iter;
    CPU_story_tot(i) = T_ADMM;

%     fprintf('Total ADMM iter: %i\n', iter)
%     fprintf('ADMM res: %i\n', res)
%     fprintf('CPU time: %f\n', T_ADMM)
end
fprintf('Average ADMM iter: %i\n', round(sum(num_iter_story_tot) / num_initial_guess))
fprintf('Average CPU time: %f\n', sum(CPU_story_tot) / num_initial_guess)

fprintf('\n')
% aAR-FP
fprintf('Starting aAR-FP.\n')

% norm_res_story_tot = zeros(maxit, num_initial_guess);
% num_iter_story_tot = zeros(num_initial_guess, 1);
% CPU_story_tot = zeros(num_initial_guess, 1);

for i = 1 : 5
    fp_iter = FPit_values(i);
    aa_iter = AAit_values(i);
    w_aAR = window_values(i);
    
    for n_guess = 1 : num_initial_guess
%         y0i = y0;
        y0i = y0(:, n_guess);
%         lambda0i = lambda0;
        lambda0i = lambda0(:, n_guess);

        x0i = [y0i; lambda0i];

        ADMMfun = @(z, u)ADMM(z, u, tol, 1, x_update, y_update, lambda_update);

        fpiterfun = @(x)fixed_point(x, ADMMfun);

        tic
        [x, iter_aAR_FP, relres_aAR, norm_story_aAR_FP] = aAR_FP(x0i, fp_iter, aa_iter, w_aAR, tol, maxit, fpiterfun);
        T_aAR_FP=toc;
        norm_res_story_tot(1:iter_aAR_FP, n_guess) = norm_story_aAR_FP;
        num_iter_story_tot(n_guess) = iter_aAR_FP;
        CPU_story_tot(n_guess) = T_aAR_FP;

%     fprintf('Total aAR-FP iter: %i\n', iter_aAR_FP)
%     fprintf('AA res: %i\n', relres_aAR)
%     fprintf('CPU time: %f\n', T_aAR_FP)
%     fprintf('\n')
    end
fprintf('aAA(%i)[%i]-FP[%i]\n',w_aAR,aa_iter,fp_iter)
fprintf('Average aAR-FP iter: %i\n', round(sum(num_iter_story_tot) / num_initial_guess))
fprintf('Average CPU time: %f\n', sum(CPU_story_tot) / num_initial_guess)
fprintf('\n')
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
