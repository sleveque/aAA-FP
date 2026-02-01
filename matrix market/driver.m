global P_l

list_string_A = ["fidap008.mtx", "fidap029.mtx", "fidapm37.mtx"];
list_string_b = ["fidap008_rhs1.mtx", "fidap029_rhs1.mtx", "fidapm37_rhs1.mtx"];

maxit_tot = 5000;
tol = 1.e-8;

fp_iter = 1;
aa_iter = 0;
w_aAR = 5000;
restart = w_aAR;
max_it = maxit_tot / restart;

fprintf("fp_iter: %i\n",fp_iter)
fprintf("a_iter: %i\n",aa_iter)
fprintf("window AA: %i\n",w_aAR)
fprintf("restart: %i\n",restart)
fprintf("maxit_tot: %i\n",maxit_tot)

for i = 1 : 3
    string_A = list_string_A(i);
    string_b = list_string_b(i);

    fprintf("\n")
    fprintf("\n")
    fprintf("problem: %s\n", string_A)

    A = mmread(string_A);
    b = mmread(string_b);

    size_A = size(A,1);
    cond_A = cond(full(A));
    cond_A_prec = cond(full(diag(A).\A));
    symm = max(max(abs(A-A')));
    fprintf("dimension of A: %i\n",size_A)
    fprintf("condition number of A: %e\n",cond_A)
    fprintf("condition number of D^{-1}A: %e\n",cond_A_prec)
    fprintf("symmetric: %e\n",full(symm))

    P_l = tril(A);

    n = size(A,1);

    x0 = ones(size(b));

    fpiterfun = @(x)fpiter(x, A, b);

    fprintf('\n')
    % GMRES
    fprintf('Starting GMRES.\n')
    tic
%    [x_GMRES, niter, relres] = my_fgmres(diag(A) .\ A, diag(A) .\ b, tol, max_it, restart, [], x0);
    [x_GMRES, niter, relres] = my_fgmres(P_l \ A, P_l \ b, tol, max_it, restart, [], x0);
    T_GMRES = toc;
    fprintf('Relative residual for GMRES: %e.\n', relres)
    fprintf('Number of GMRES iterations: %i.\n', (niter(1,1)-1)*restart+niter(1,2))
    fprintf('CPU time: %f\n', T_GMRES)

    fprintf('\n')
    % AA
    fprintf('Starting AA.\n')
    tic
    [x_aAR_FP, iter_AA, relres_AA, anorm_story_aAR_FP, rnorm_sotry_aAR_FP] = aAR_FP(x0, 0, 1, w_aAR, tol, maxit_tot, fpiterfun);
    T_AA = toc;
    fprintf('Total AA iter: %i\n', iter_AA)
    fprintf('AA res: %i\n', relres_AA)
    fprintf('CPU time: %f\n', T_AA)

    fprintf('\n')
    % aAR-FP
    fprintf('Starting aAR-FP.\n')
    tic
    [x, iter_aAR_FP, relres_aAR, norm_story_aAR_FP] = aAR_FP(x0, fp_iter, aa_iter, w_aAR, tol, maxit_tot, fpiterfun);
    T_aAR_FP=toc;

    fprintf('Total aAR-FP iter: %i\n', iter_aAR_FP)
    fprintf('AA res: %i\n', relres_aAR)
    fprintf('CPU time: %f\n', T_aAR_FP)

end


function xnew = fpiter(xold, A, b)
    global P_l

    xnew = b - A * xold;
%    omega = 0.5;
    
    xnew = xold + P_l \ xnew;
%    xnew = xold + omega * diag(A) .\ xnew;
end
