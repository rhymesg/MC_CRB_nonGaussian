# Numerical regression checks

Closed-form Gaussian and Laplace values test central and far-tail log densities, including density-underflow cases.

From the repository root, with base MATLAB:

```bash
matlab -batch "addpath('tests/integration/likelihood'); verify_likelihood"
```

These checks require no external data or plotting. They have been syntax checked, but have not been executed in MATLAB or Octave. They do not validate the complete research experiment.
