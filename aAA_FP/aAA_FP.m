function [xnew, iter, tol, anorm_story, rnorm_story, x_story] = aAA_FP(x0, it_fp, it_aa, m, rtol, maxit, fpiter)

% aAA(m)[it_aa]-FP[it_fp] iteration
%
% input
%
% x0             intial guess
% it_fp          number of fixed point iterations
% it_aa          number of Anderson acceleration
% m              window of AA
% rtol           relative tolerance on the residual to achieve
% maxit          number of maximum iterations
% fpiter         fixed point iteration to accelerate
%
% output
%
% xnew           the solution evaluated
% iter           number of total iterations required
% tol            tolerance achieved
% anorm_story    story of the absolute tolerances
% rnorm_story    story of the relative tolerances
% x_story        story of the solutions at each iterate


    % vectors containing the story of the absolute and relative tolerances
    anorm_story = zeros(maxit+1,1);
    rnorm_story = zeros(maxit+1,1);

    % the first iteration is always a fixed point
    xnew = fpiter(x0);
    rold = xnew - x0;
    norm0 = norm(rold);

    anorm_story(1) = norm0;
    rnorm_story(1) = 1;

    % vector containing
    x_story(:,1) = x0;

    % matrices containing the windowed story of residuals and solutions
    Rrk = rold;
    Qk = xnew;
    
    for k = 2 : maxit
        x_story(:,k) = xnew;

        % applying the fixed point iteration
        xfp = fpiter(xnew);
        rnew = xfp - xnew;
        
        norm_k = norm(rnew);
        
        anorm_story(k) = norm_k;
        rnorm_story(k) = norm_k / norm0;

        if norm_k / norm0 < rtol
            break
        end

        % checking the dimension of Rrk and Qk
        mk = min(m,k-1);

        % updating Rrk and Qk
        e = ones(1,mk);
        Rrk = [rnew Rrk(:,1:mk)];
        Qk = [xfp Qk(:,1:mk)];
        Rk = Rrk(:,2:end)-kron(e,rnew);
        Xk = kron(e,xfp) - Qk(:,2:end);

        % test to understand whether employing FP or applying AA
        mod_k = mod(k - 1, it_fp + it_aa);
        
        if mod_k < it_fp
            % fixed point iteration
            fprintf('iter: %i, fixed point\n',k)
            xnew = xfp;
        else
            % Anderson acceleration
            fprintf('iter: %i, Anderson Acceleration\n',k)
            coeff = Rk \ rnew;

            xnew = xfp + Xk * coeff;
        end

    end

    iter = k;
    tol = norm(rnew) / norm0;
    anorm_story=anorm_story(1:k);
    rnorm_story=rnorm_story(1:k);

end
