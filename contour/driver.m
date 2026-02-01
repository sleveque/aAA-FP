global dim_y dim_lambda x_sol y_sol lambda_sol

seed = 100000;
rng(seed)

tol = 1e-12;
maxit = 1000;


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

ADMMfun = @(z, u)ADMM(z, u, tol, 1, x_update, y_update, lambda_update);

fpiterfun = @(x)fixed_point(x, ADMMfun);

CPU = zeros(20, 20);
it = zeros(20, 20);

for n_range = 1 : 20
    for m_range = 1 : 20
        fprintf('\n')
        % aAR-FP
        fprintf('Starting aAR-FP.\n')

        % window for AA
        w_aAR = n_range;
        % number of AA iterations
        aa_iter = n_range;
        % number of FP iterations
        fp_iter = m_range;

        fprintf('AA window: %i\n', w_aAR)
        fprintf('Number of AA iterations: %i\n', aa_iter)
        fprintf('Number of FP iterations: %i\n', fp_iter)

        x0 = [y0; lambda0];

        tic
        [x, iter_aAR_FP, relres_aAR, norm_story_aAR_FP] = aAR_FP(x0, fp_iter, aa_iter, w_aAR, tol, maxit, fpiterfun);
        T_aAR_FP=toc;

        fprintf('Total aAR-FP iter: %i\n', iter_aAR_FP)
        fprintf('AA res: %i\n', relres_aAR)
        fprintf('CPU time: %f\n', T_aAR_FP)
        CPU(n_range, m_range) = T_aAR_FP;
        it(n_range, m_range) = iter_aAR_FP;
    end
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
