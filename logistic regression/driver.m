global y v mu n eta

addpath('libsvm-3.36/matlab')

eta = 1;
mu = 1e-02;

tol = 1e-12;
maxit = 1000;

m_AA = 5;
it_aa = 5;
it_fp = 10;

dim_x = 54;
n = 581012;

[y, v] = libsvmread('covtype.libsvm.binary.scale');
v = v';

% gradient descent
xold = zeros(dim_x, 1);
res_GD = zeros(maxit, 1);

tic
for iter = 1 : maxit
    xnew = fp_iter(xold);
    res = norm(xold - xnew);
    res_GD(iter) = res;
    fprintf('iter %i, norm: %e\n', iter, res)
    if res < tol
        break
    end
    xold = xnew;
end
time_GD = toc;

fprintf('Total GD iter: %i\n', iter)
fprintf('GD res: %i\n', res)
fprintf('CPU time: %f\n', time_GD)


% aAA(m)[s]-FP[t]
xold = zeros(dim_x, 1);

fp_iter_fun = @(x)fp_iter(x);
tic
[xnew, iter, tol, anorm_story, rnorm_story, x_story] = aAR_FP(xold, it_fp, it_aa, m_AA, tol, maxit, fp_iter_fun);
time_aAA_FP = toc;

fprintf('Total aAA-FP iter: %i\n', iter)
fprintf('aAA-FP res: %i\n', tol)
fprintf('CPU time: %f\n', time_aAA_FP)


function xnew = fp_iter(x)
    global eta

    xnew = x - eta * gradient_eval(x);
end


function grad_h = gradient_eval(x)
    global y v mu n

    M_help = y' .* v;
    v_help = v' * x;
    v_help = exp(- v_help .* y);
    v_help = v_help ./ (1 + v_help);
    grad_h = - M_help * v_help;
    grad_h = (1/n) * grad_h;
    grad_h = grad_h + mu * x;

end
