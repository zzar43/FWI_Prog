# Frequency Domain Full Waveform Inversion

This is a new version frequency domain full waveform inversion (FWI) code in [Julia](https://julialang.org).

A 5 point finite difference solver is used to solve the Helmholtz equation with perfect matched layer.
And the l-bfgs-b method is used during the optimization which is provided by https://github.com/Gnimuc/LBFGSB.jl.

For demo please check [Reduced_Prob.ipynb](https://github.com/zzar43/FWI_Prog/blob/master/Reduced_Prob.ipynb).

Next:
- An iterative helmholtz solver with inexact LU decomposition used as the preconditioner.
- Distributed computing during the multi-sources forward modelling.
- Add more regularization as TV, develop new regularization techniques.
- New optimization technique such as penalty method, consider application on the multi-variable FWI.