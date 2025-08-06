function [xnew, iter, tol, anorm_story, rnorm_story] = aAA_FP(x0, it_fp, it_aa, m, rtol, maxit, fpiter)

    anorm_story = zeros(maxit+1,1);
    rnorm_story = zeros(maxit+1,1);
    
    xnew = fpiter(x0);
    rold = xnew - x0;
    norm0 = norm(rold);

    anorm_story(1) = norm0;
    rnorm_story(1) = 1;
    
    Rrk = rold;
    Qk = xnew;
    
    for k = 2 : maxit
        
        xfp = fpiter(xnew);
        rnew = xfp - xnew;
        
        norm_k = norm(rnew);
        
        anorm_story(k) = norm_k;
        rnorm_story(k) = norm_k / norm0;

        if norm_k < rtol
            break
        end

        mk = min(m,k-1);

        e = ones(1,mk);
        Rrk = [rnew Rrk(:,1:mk)];
        Qk = [xfp Qk(:,1:mk)];
        Rk = Rrk(:,2:end)-kron(e,rnew);
        Xk = kron(e,xfp) - Qk(:,2:end);
        
        mod_k = mod(k - 1, it_fp + it_aa);
        
        if mod_k < it_fp
             fprintf('iter: %i, fixed point\n',k)
            xnew = xfp;
        else
             fprintf('iter: %i, Anderson Acceleration\n',k)
            coeff = Rk \ rnew;

            xnew = xfp + Xk * coeff;
        end

    end

    iter = k;
    tol = norm(rnew);
    anorm_story=anorm_story(1:k);
    rnorm_story=rnorm_story(1:k);

end
