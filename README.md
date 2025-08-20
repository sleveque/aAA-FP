Maltab code for the numerical examples in the following paper:
[Add arXiv link](https://arxiv.org/abs/2508.10158)

@article{aAAFP25,
	title={A Generalized Alternating Anderson Acceleration Method},
	author={He, Yunhui and Leveque, Santolo},
	journal={arXiv preprint arXiv:2508.10158},
	year={2025}
}

The code aAA_FP implements the generalized alternating Anderson acceleration method aAA($m$)[ $s$ ]–FP[ $t$ ] to accelerate fixed-point iteration, which is a periodic
scheme composed of $t$ fixed-point iteration steps, interleaved with $s$ steps of Anderson acceleration
with window size m, to solve linear and nonlinear problems. Specifically, given the problem

$g(x)=0$,

whose solution $x^*$ solves the fixed-point iteration

$q({x^*}) = {x^ *}$,

aAA_FP applies $s$ iterations of AA($m$) after performing $t$ fixed-point iterations,
then repeat this process until a reduction on the residual is achieved.

The call to the routine is the following:

[xnew, iter, tol, anorm_story, rnorm_story, x_story] = aAA_FP(x0, it_fp, it_aa, m, rtol, maxit, fpiter)


Input:

- x0             intial guess
- it_fp          number of fixed point iterations
- it_aa          number of Anderson acceleration
- m              window of AA
- rtol           relative tolerance on the residual to achieve
- maxit          number of maximum iterations
- fpiter         fixed point iteration to accelerate

Output:

- xnew           the solution evaluated
- iter           number of total iterations required
- tol            tolerance achieved
- anorm_story    story of the absolute tolerances
- rnorm_story    story of the relative tolerances
- x_story        story of the solutions at each iterate

Tests: we apply aAA-FP to accelerate several types of fixed-point iterations:       

        1. Jacobi iteration for $solving Ax=b$, where the coefficient matrices are from the Matrix Market repository.
        2. Picard iteration for Navier-Stokes equations. This code is available from the authors upon request.
        3. Alternating direction method of mulitpliers method (ADMM).
        4. Gradient descent method for regularized logistrtic regression.
        
