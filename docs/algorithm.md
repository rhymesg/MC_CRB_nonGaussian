# Monte Carlo information and terrain-navigation bounds

This reference describes the MATLAB implementation accompanying the [APISAT 2017 paper](../README.md#citation), based on source revision [`fd7df2f`](https://github.com/rhymesg/MC_CRB_nonGaussian/tree/fd7df2f1899965aee68adc1c58ce78a3d5e8d77e). It supports reading and adapting the code; equation-level agreement with the unavailable paper text has not been established.

## Measurement and likelihood contracts

The simulated observation is `z = h(x) + noise`, with synthetic terrain `h(x) = 70 + 12 sin(2πx/80)`. Position and height have no declared physical unit in the source; retain consistent position, height, and noise units when adapting it.

| Function | Inputs | Output and assumptions |
|---|---|---|
| [meas_TRN_1d](../meas_TRN_1d.m) | Scalar or array `x` | Same-shaped terrain height; sinusoid is defined beyond the plotted interval |
| [loglikelihood_TRN](../loglikelihood_TRN.m) | Scalar predicted/observed heights `Z_est`, `Z`, positive standard deviation `sig_Z`, selector `Dist` | Scalar natural log density; `Dist == 0` is Gaussian, every other value selects Laplace |
| [laprnd](../laprnd.m) | Dimensions `m,n`; optional mean `mu` and standard deviation `sigma` | `m × n` samples; defaults are zero mean and unit standard deviation; Laplace scale is `sigma/sqrt(2)` |

For residual `r = Z - Z_est` and `σ = sig_Z`, the intended log densities are:

- Gaussian: `−log(σ√(2π)) − r²/(2σ²)`.
- Laplace: `−log(2b) − |r|/b`, where `b = σ/√2`.

The helper computes the density first and then takes `log`; extreme residuals can underflow to `-Inf`. Inputs are not comprehensively validated, and the scalar likelihood's matrix-power operators must not be treated as a vectorized API.

## Simulation entry points

- [main_TRN_1d.m](../main_TRN_1d.m) redraws particles at each position with Gaussian offsets and weights them using the observation likelihood; it does not propagate a posterior between positions.
- [main_TRN_1d_recur.m](../main_TRN_1d_recur.m) propagates particles by the known trajectory increment, a shared Gaussian perturbation, and independent Gaussian perturbations; it accumulates likelihood weights and resamples when effective sample size falls below `0.3*Np`.
- `sig_Z` controls measurement standard deviation; `sig` and `sig_p` control Gaussian state offsets and process perturbations. `Np` sets particles, `N_MC` sets filter trials, `L` sets sampled states, and `N` sets pseudodata per sampled state.
- `x_est` and `x_err` are `N_MC × LENG`; RMSE and bound outputs are row vectors of length `LENG`. Recursive `FIM_*_res(1)` is zero from MATLAB array growth, not an estimated first-step measurement information value.

## Information estimation

Both scripts compare these calculations for sampled states near each trajectory position:

1. Approximate terrain slope with a centered difference: `H = (h(θ+s) − h(θ−s))/(2s)`.
2. Average `H²/sig_Z²` across states for `FIM_lin`. This is a Gaussian-form comparison even when observations use Laplace noise.
3. For each state `θ`, draw pseudodata `z` from its selected observation distribution and independent signs `Δ, Δ̃ ∈ {−1,+1}`.
4. Evaluate `g± = [ℓ(θ ± cΔ + c̃Δ̃; z) − ℓ(θ ± cΔ; z)]/(c̃Δ̃)`, where `ℓ(θ;z)` calls the terrain model and log-likelihood helper.
5. Estimate the scalar Hessian as `(g+ − g−)/(2cΔ)`, average over pseudodata, and negate the result.
6. If the result is negative, the code applies `sqrtm(F*F)` before averaging over states; for a real scalar this is `abs(F)`.

In the independent-position script, `J = 1/(sig²+sig_p²) + FIM`. In the recursive script, each comparison starts at `J = 1/sig²` and updates `J = 1/(1/J + sig_p²) + FIM`; the reported bound is `sqrt(1/J)`.

## Numerical interpretation and adaptation

- Recursive state sampling uses `prior = inv(J)*randn(1,L)`: the multiplier is inverse information, rather than its square root. Preserve this distinction when translating; changing it changes the experiment.
- Laplace log density has a cusp at zero residual; finite-difference Hessian estimates depend on perturbation sizes and sampled residuals. Absolute-value correction of negative estimates is an implementation choice, not evidence of an unbiased FIM estimate.
- Particle weights are products of densities without log-space normalization; zero total weight can produce `NaN` values.
- Internal `rng('shuffle')` calls prevent a caller-supplied seed from fixing the full run. A repeatable adaptation needs controlled random draws throughout, and matching seeds across languages does not imply matching samples.
- Preserve MATLAB's one-based time indexing, row-vector shapes, and multiplication order when translating to Python or C++. The scalar sign/absolute-value operation should not be generalized into a matrix positive-semidefinite projection without a separate derivation.
- Compare a translation first against the [deterministic helper example](../README.md#usage), then against intermediate likelihood and Hessian calculations using identical supplied random draws. Full simulation curves have no recorded reference values or justified comparison tolerances here.

The [related toolkit](https://github.com/rhymesg/information-based-tracking) uses score outer products for Monte Carlo information. That alternative and these approximate recursive updates must not be assumed interchangeable or established as an exact Bayesian bound for arbitrary non-Gaussian models.
