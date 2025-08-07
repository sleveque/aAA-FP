The code aAA_FP implements the generalized alternating Anderson acceleration method aAA($m$)[ $s$ ]–FP[ $t$ ]. Specifically, given the problem

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
