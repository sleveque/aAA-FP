list_string_A = ["fidap008.mtx", "fidap029.mtx", "fidapm37.mtx"];
list_string_b = ["fidap008_rhs1.mtx", "fidap029_rhs1.mtx", "fidapm37_rhs1.mtx"];

maxit_tot = 10000;
tol = 1.e-8;

fp_iter = 4;
aa_iter = 1;
w_aAA = 10;
restart = 100;
max_it = maxit_tot / restart;

fprintf("fp_iter: %i\n",fp_iter)
fprintf("a_iter: %i\n",aa_iter)
fprintf("window AA: %i\n",w_aAA)
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

    n = size(A,1);

    x0 = ones(size(b));

    fpiterfun = @(x)fpiter(x, A, b);

    mfun1 = @(x)prec(x, A, b);

    fprintf('\n')
    % AA
    fprintf('Starting AA.\n')
    tic
    [x_aAA_FP, iter_AA, relres_AA, anorm_story_aAA_FP, rnorm_sotry_aAA_FP] = aAA_FP(x0, 0, aa_iter, w_aAA, tol, max_it * restart, fpiterfun);
    T_AA = toc;
    fprintf('Total AA iter: %i\n', iter_AA)
    fprintf('AA res: %i\n', relres_AA)
    fprintf('CPU time: %f\n', T_AA)

    fprintf('\n')
    % GMRES
    fprintf('Starting GMRES.\n')
    tic
    % restarted FGMRES
    % [x_GMRES, niter, relres] = my_fgmres(A, b, tol, max_it, restart, mfun1, x0);
    % full FGMRES (no restart)
    [x_GMRES, niter, relres] = my_fgmres(diag(A).\A, diag(A).\b, tol, max_it, restart, [], x0);
    T_GMRES = toc;
    fprintf('Relative residual for GMRES: %e.\n', relres)
    fprintf('Number of GMRES iterations: %i.\n', (niter(1,1)-1)*restart+niter(1,2))
    fprintf('CPU time: %f\n', T_GMRES)

    fprintf('\n')
    % aAA_FP
    fprintf('Starting aAA_FP.\n')
    tic
    [x, iter_aAA_FP, relres_aAA, norm_story_aAA_FP] = aAA_FP(x0, fp_iter, aa_iter, w_aAA, tol, max_it * restart, fpiterfun);
    T_aAA_FP=toc;

    fprintf('Total aAA-FP iter: %i\n', iter_aAA_FP)
    fprintf('AA res: %i\n', relres_aAA)
    fprintf('CPU time: %f\n', T_aAA_FP)

end


function xnew = fpiter(xold, A, b)
    xnew = b - A * xold;
    
    xnew = diag(A) .\ xnew;
end

function xnew = prec(xold, A, b)
    D = diag(A);

    res = b - A * xold;

    xnew = D .\ res;
end
