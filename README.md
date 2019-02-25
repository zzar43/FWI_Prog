# Frequency Domain Full Waveform Inversion

This is a new version frequency domain full waveform inversion (FWI) code.
A 5 point finite difference solver is used to solve the Helmholtz equation with perfect matched layer.
And the l-bfgs-b method is used during the optimization which is provided by https://github.com/Gnimuc/LBFGSB.jl.

For demo please check Reduced_Prob.ipynb.