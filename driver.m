global y v mu n eta

addpath('libsvm-3.36/matlab')
num_test = 2;

eta = 1;
mu = 1e-02;

tol = 1e-12;
maxit = 100;

m_AA = 3;
it_aa = 1;
it_fp = 0;


if num_test == 1
    dim_x = 300;
    n = 49749;

    [y, v] = libsvmread('w8a');
    v = v';
end

if num_test == 2
    dim_x = 54;
    n = 581012;

    [y, v] = libsvmread('covtype.libsvm.binary');
    v = v';
end

% % gradient descent
% xold = zeros(dim_x, 1);
% 
% tic
% for iter = 1 : maxit
%     xnew = fp_iter(xold);
%     res = norm(xold - xnew);
%     fprintf('iter %i, norm: %e\n', iter, res)
%     if res < tol
%         break
%     end
%     xold = xnew;
% end
% time_GD = toc;
% 
% fprintf('Total GD iter: %i\n', iter)
% fprintf('GD res: %i\n', res)
% fprintf('CPU time: %f\n', time_GD)

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

    grad_h = zeros(size(x));

    dim_x = size(x, 1);

    for j = 1 : dim_x
        for i = 1 : n
            d_dxj = - y(i) * v(j, i) * exp(- y(i) * (x' * v(:, i)));
            d_dxj = d_dxj / (1 + exp(- y(i) * (x' * v(:, i))));
            grad_h(j) = grad_h(j) + d_dxj;
        end
        grad_h(j) = (1 / n) * grad_h(j);
        grad_h(j) = grad_h(j) + mu * x(j);
    end
end
