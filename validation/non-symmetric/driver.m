maxit = 100;
tol = 1.e-10;

fp_iter = 0;
aa_iter = 1;
w_aAA = maxit;

n = 32;

b = zeros(n,1);
b(1,1) = 1;

A = [zeros(1,n-1), 1; eye(n-1), zeros(n-1,1)];

x0 = ones(size(b));

% aAA_FP
fpiterfun = @(x)fpiter(x, A, b);

[x_aAA_FP, iter_aAA_FP, relres_aAA, anorm_story_aAA_FP, rnorm_sotry_aAA_FP, x_story_aAA_FP] = aAA_FP(x0, fp_iter, aa_iter, w_aAA, tol, maxit, fpiterfun);

% story of solutions of GMRES
x_GMRES = [x0];
for n_iter = 1 : n
    [x_GMRES_i, e, res] = my_gmres(A, b, x0, n_iter, tol);
    
    x_GMRES = [x_GMRES, x_GMRES_i];
    
    clear x_GMRES_i
end

% evaluating residuals of GMRES
res_story_GMRES = [];
for i = 1 : n + 1
    res_story_GMRES = [res_story_GMRES, reseval(x_GMRES(:, i), A, b)];
end

% applying M to GMRES residuals
fp_res_GMRES = [];
for i = 1 : n + 1
    fp_res_GMRES = [fp_res_GMRES, (eye(size(A)) - A) * res_story_GMRES(:, i)];
end

% evaluating residuals of aAA_FP
res_story_aAA_FP = [];
for i = 1 : iter_aAA_FP
    res_story_aAA_FP(:, i) = reseval(x_story_aAA_FP(:, i), A, b);
end

% evaluating norms
norm_fpiter_GMRES = zeros(n,1);
for i = 1 : n + 1
    norm_fpiter_GMRES(i, 1) = norm(fp_res_GMRES(:,i));
end

norm_res_aAA_FP = zeros(iter_aAA_FP, 1);
for i = 1 : iter_aAA_FP
    norm_res_aAA_FP(i, 1) = norm(res_story_aAA_FP(:, i));
end

% plots
figure
ax=gca;
ax.FontSize = 14;
set(gca,'TickLabelInterpreter','latex')
axis([-1 length(norm_res_aAA_FP) min(min(norm_fpiter_GMRES),min(norm_res_aAA_FP)) max(norm_res_aAA_FP)])
xlabel('$k$','interpreter','latex','FontSize',18)
ylabel('$\|r_{k}\|$','interpreter','latex','FontSize',18)
hold

x=linspace(0, length(norm_fpiter_GMRES)-1, length(norm_fpiter_GMRES));
semilogy(x,norm_fpiter_GMRES, 'x-')
x=linspace(-1, length(norm_res_aAA_FP)-2, length(norm_res_aAA_FP));
semilogy(x,norm_res_aAA_FP, 'o-')
legend('$\|M r_{k}^{G}\|$', '$\|r_{k}\|$','interpreter','latex','FontSize',18)



function res = reseval(x, A, b)
    res = b - A * x;
end

function xnew = fpiter(xold, A, b)
    res = reseval(xold, A, b);
    
    xnew = xold + res;
end
