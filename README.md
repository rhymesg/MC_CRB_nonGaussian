# MC_CRB_nonGaussian

## Overview

MATLAB simulations of Monte Carlo Fisher information matrix (FIM) estimation and Cramér–Rao bounds (CRB) for one-dimensional terrain-referenced navigation with Gaussian or Laplace measurement noise.

The [canonical repository](https://github.com/rhymesg/MC_CRB_nonGaussian) accompanies the [APISAT 2017 paper](#citation). Compare a linearized measurement-information calculation with a simultaneous-perturbation, log-likelihood Hessian estimate and particle-filter RMSE.

For reusable measurement models, information prediction, and a **score-based** Monte Carlo estimator, see [information-based-tracking](https://github.com/rhymesg/information-based-tracking). Its [estimator](https://github.com/rhymesg/information-based-tracking/blob/main/monte_carlo_information.m) uses score outer products, providing a complementary information-estimation approach.

For Python, C++, or other-language implementations, use the [algorithm and translation reference](docs/algorithm.md) to follow the likelihood formulas, array contracts, and information-update order.

## Method

Estimate measurement information in two ways: a local linearized model and a Monte Carlo finite-difference Hessian of the log likelihood. The scripts combine these calculations with a one-dimensional terrain particle filter to compare uncertainty bounds and estimation error.

### Implementation reference

| Purpose | Source |
|---|---|
| Recursive particle filter and information prediction/update | [main_TRN_1d_recur.m](main_TRN_1d_recur.m) |
| Independent-position particle estimates and bound comparison | [main_TRN_1d.m](main_TRN_1d.m) |
| Synthetic sinusoidal terrain measurement | [meas_TRN_1d.m](meas_TRN_1d.m) |
| Scalar Gaussian/Laplace log likelihood | [loglikelihood_TRN.m](loglikelihood_TRN.m) |
| Laplace sampling with a standard-deviation parameter | [laprnd.m](laprnd.m) |

The [algorithm reference](docs/algorithm.md) documents the MATLAB inputs, update order, and numerical conventions.

## Examples

Run from the repository root in MATLAB. The source uses base MATLAB and synthetic terrain; no additional toolbox or dataset is needed. Shell commands require MATLAB's `-batch` option.

Run the recursive particle-filter and bound example:

```bash
matlab -batch "main_TRN_1d_recur"
```

Run the independent-position comparison:

```bash
matlab -batch "main_TRN_1d"
```

Edit `Dist` inside the chosen script: `0` selects Gaussian measurement noise and `1` selects Laplace noise. The recursive script defaults to Laplace; the independent-position script defaults to Gaussian.

Each script clears the workspace, runs nested Monte Carlo loops, prints `PF completed` before estimating information, and plots terrain, particle-filter RMSE, and two bound curves. Workspace outputs include `x_err_RMS`, `FIM_lin_res`, `FIM_mon_res`, `LB1_res`, and `LB2_res`; the `LB` arrays contain standard-deviation bounds, not variances.

For a small deterministic example, run this in MATLAB from the repository folder:

```matlab
x = [0, 20, 40];
height = meas_TRN_1d(x)
gaussian_log_density = loglikelihood_TRN(70, 70, 2, 0)
laplace_log_density = loglikelihood_TRN(70, 70, 2, 1)
```

Expected values, derived from the helper formulas: `height = [70, 82, 70]`, Gaussian log density approximately `-1.6120857138`, and Laplace log density approximately `-1.0397207708`.

## Implementation scope

The helper routines expose terrain measurements, Gaussian/Laplace log likelihoods, and Laplace sampling. The simulation scripts show particle-filter and information-recursion calculations; [calculation details](docs/algorithm.md) connect their formulas, parameters, and source behavior.

### Checks

Run the focused [Gaussian/Laplace log-density checks](tests/integration/likelihood/README.md):

```bash
matlab -batch "addpath('tests/integration/likelihood'); verify_likelihood"
```

The full simulations initialize random streams with `rng('shuffle')` inside their loops; control those calls when constructing a repeatable experiment.

## Citation

Please cite the related paper when using its method:

Youngjoo Kim and Hyochoong Bang. “Monte-Carlo Calculation of Cramer-Rao Bound for non-Gaussian Recursive Filtering.” *2017 Asia-Pacific International Symposium on Aerospace Technology* (APISAT), 2017.

[CITATION.cff](CITATION.cff) provides machine-readable metadata. Publication details follow the original author README.

## License and provenance

The repository includes an [MIT license](LICENSE). Preserve its copyright and permission notices when reusing the software; paper citation is a separate academic attribution request.

The FIM calculation is based on Sonjoy Das's [estimateFIM.html](estimateFIM.html); `laprnd.m` retains Elvis Chen's original attribution. See [provenance](docs/provenance.md) for the upstream reference, citation request, and third-party licensing information.
