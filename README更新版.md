# VFCR-ZNN Example 1

Place the following files in the same MATLAB directory:

- `run_example1.m`
- `example1_problem.m`
- `my_system.m`
- `znn_ode.m`
- `vfcr_residual.m`

Run `run_example1.m`.

The ODE state is

```matlab
Y = [g; z]
```

where `g=[x;mu1;mu2]` has 7 entries and the history-integral state `z`
also has 7 entries. Therefore Equation (20) is integrated as a
14-dimensional ODE.

The default `phi_eps=1e-4` uses a smooth approximation to the NFTAF only
near zero. This prevents `ode15s` from taking extremely small steps at the
non-Lipschitz point. Set `phi_eps=0` in `run_example1.m` for the literal
Equation (18), but expect a much longer runtime.

The paper does not report its random seed, ODE solver, tolerances, or
numerical treatment at zero. The reproducible target is therefore the
qualitative fixed-time behavior and agreement with the time-varying KKT
solution, not pixel-for-pixel equality of the initial transients.
