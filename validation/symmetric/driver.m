% DEFINITION OF THE INPUT
a=-1;b=1;c=-1;d=1;
% order of finite element discretization
FE_order = 2;

% linear solver parameters
maxit = 100;
tol = 1.e-12;

% aAA-FP parameters
fp_iter = 3;
aa_iter = 1;
w_aAA = maxit;


warning off

nx_exp=4;

for boh_2 = 1 : 1
nx_exp=nx_exp+1;
ny_exp=nx_exp;

n = 2 ^ nx_exp - 1;
m = 2 ^ ny_exp - 1;

n_space = n * m;

% mesh-size in space
h = max((b - a) / (n + 1),(d - c) / (m + 1));

fprintf('\n')
fprintf('grid refinement level: %i\n', nx_exp)
fprintf('order FE: %i\n', FE_order)
fprintf('number of points on the x-axis: 2^%i+1\n', nx_exp)
fprintf('mesh-size in space: %f\n', h)

% the code can create finite element matrices also for Shishkin (rectangular) meshes
sigma_x=0;
sigma_y=0;

number_nodes=2;
[K, K_boundary, rhs] = FEM_q1_discretize(a, b, c, d, nx_exp, ny_exp, number_nodes, sigma_x, sigma_y);

v_boundary_K = construct_boundary_conditions_for_v(a, b, c, d, nx_exp, ny_exp, sigma_x, sigma_y, K_boundary, FE_order);

rhs = rhs - v_boundary_K;

fprintf('dimension of the system solved: %i\n', n_space)

true_xsol = zeros(n_space,1);

if sigma_x == 0 && sigma_y == 0
    x1=linspace(a,b,n+2);
    x2=linspace(c,d,m+2);
else
    x1 = zeros(n+2,1);
    x2 = zeros(m+2,1);

    x1(1:(n+1)/2+1) = linspace(a,sigma_x,(n+3)/2);
    x1((n+1)/2+1:end) = linspace(sigma_x,b,(n+3)/2); % overwrite middle point
    x2(1:(m+1)/2+1) = linspace(c,sigma_y,(m+3)/2);
    x2((m+1)/2+1:end) = linspace(sigma_y,d,(m+3)/2); % overwrite middle point
end

x1_1=x1(2:n+1);
x2_1=x2(2:m+1);

for i=1:n
    for j=1:m
        true_xsol((i - 1) * m + j, 1) = cos(pi*x1_1(i)/2)*cos(pi*x2_1(j)/2) + 1;
    end
end

x0 = zeros(size(rhs));

% exact solver
xsol = K \ rhs;

% aAA_FP
fpiterfun = @(x)fpiter(x, K, rhs);

[x_aAA_FP, iter_aAA_FP, relres_aAA, anorm_story_aAA_FP, rnorm_sotry_aAA_FP, x_story_aAA_FP] = aAA_FP(x0, fp_iter, aa_iter, w_aAA, tol, maxit, fpiterfun);

% GMRES
x_GMRES = zeros(size(rhs,1), maxit + 1);
x_GMRES(:, 1) = x0;
for n_iter = 1 : maxit
    [x_GMRES_i, e, res] = my_gmres(K, rhs, x0, n_iter, tol);
    
    x_GMRES(:, n_iter + 1) = x_GMRES_i;

    if e(end) <= tol
        break
    end
end

it = size(e, 1);
x_GMRES = x_GMRES(:, 1 : it);

% evaluating residual of GMRES
res_story_GMRES = [];
for i = 1 : it
    res_story_GMRES = [res_story_GMRES, reseval(x_GMRES(:, i), K, rhs)];
end

% applying M
fp_res_GMRES = [];
for i = 1 : it
    fp_res_GMRES = [fp_res_GMRES, (eye(size(K)) - K) * res_story_GMRES(:, i)];
end

% evaluating residuals of aAA_FP
res_story_aAA_FP = [];
for i = 1 : iter_aAA_FP
    res_story_aAA_FP(:, i) = reseval(x_story_aAA_FP(:, i), K, rhs);
end

% evaluating norms
norm_fpiter_GMRES = zeros(it,1);
for i = 1 : it
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

fprintf('\n')

end



function res = reseval(x, A, b)
    res = b - A * x;
end

function xnew = fpiter(xold, A, b)
    res = reseval(xold, A, b);
    
    xnew = xold + res;
end
