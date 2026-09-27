# MC_CRB_nonGaussian

## Overview

MATLAB simulations of Monte Carlo Fisher information matrix (FIM) estimation and Cramér–Rao bounds (CRB) for one-dimensional terrain-referenced navigation with Gaussian or Laplace measurement noise.

The [canonical repository](https://github.com/rhymesg/MC_CRB_nonGaussian) accompanies the [APISAT 2017 paper](#citation). Compare a linearized measurement-information calculation with a simultaneous-perturbation, log-likelihood Hessian estimate and particle-filter RMSE.

For reusable measurement models, information prediction, and a **score-based** Monte Carlo estimator, see [information-based-tracking](https://github.com/rhymesg/information-based-tracking). Its [estimator](https://github.com/rhymesg/information-based-tracking/blob/main/monte_carlo_information.m) uses score outer products; it does not reproduce this repository's perturbation estimator or recursive navigation simulation.

## Installation

Clone the repository:

```bash
git clone https://github.com/rhymesg/MC_CRB_nonGaussian.git
```

Enter its directory:

```bash
cd MC_CRB_nonGaussian
```

Use MATLAB with this directory as the current folder. The scripts use base MATLAB functions and synthetic terrain; no dataset or additional toolbox is required by the source. No minimum MATLAB release or Octave compatibility has been established; shell commands below require MATLAB's `-batch` option.

## Usage

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

## Development

There is no automated test suite or CI configuration. The deterministic example above checks helper behavior, not the paper's simulation results; native MATLAB and Octave execution remain unverified.

Report problems through [GitHub Issues](https://github.com/rhymesg/MC_CRB_nonGaussian/issues), including the script, parameter changes, MATLAB version, and error or unexpected output. The scripts call `rng('shuffle')` internally, so setting a seed before running them does not make the full simulation repeatable.

## Implementation reference

| Purpose | Source |
|---|---|
| Recursive particle filter and information prediction/update | [main_TRN_1d_recur.m](main_TRN_1d_recur.m) |
| Independent-position particle estimates and bound comparison | [main_TRN_1d.m](main_TRN_1d.m) |
| Synthetic sinusoidal terrain measurement | [meas_TRN_1d.m](meas_TRN_1d.m) |
| Scalar Gaussian/Laplace log likelihood | [loglikelihood_TRN.m](loglikelihood_TRN.m) |
| Laplace sampling with a standard-deviation parameter | [laprnd.m](laprnd.m) |

The [algorithm reference](docs/algorithm.md) documents inputs, update order, numerical limitations, and translation considerations for Python or C++ readers. This repository supplies MATLAB source, without ports or language bindings.

## Citation

Please cite the related paper when using its method:

Youngjoo Kim and Hyochoong Bang. “Monte-Carlo Calculation of Cramer-Rao Bound for non-Gaussian Recursive Filtering.” *2017 Asia-Pacific International Symposium on Aerospace Technology* (APISAT), 2017.

[CITATION.cff](CITATION.cff) provides machine-readable metadata. Publication details are retained from the original author README; no paper DOI or full-text link has been verified.

## License and provenance

The repository includes an [MIT license](LICENSE). Preserve its copyright and permission notices when reusing the software; paper citation is a separate academic attribution request.

The FIM calculation is based on Sonjoy Das's [estimateFIM.html](estimateFIM.html); `laprnd.m` retains Elvis Chen's original attribution. See [provenance](docs/provenance.md) for the upstream reference, citation request, and third-party licensing information.
