# Source provenance and reuse

This repository contains Youngjoo Kim's MATLAB terrain-navigation simulations associated with the [APISAT 2017 paper](../README.md#citation). The historical baseline is [`fd7df2f`](https://github.com/rhymesg/MC_CRB_nonGaussian/tree/fd7df2f1899965aee68adc1c58ce78a3d5e8d77e); current log densities are evaluated directly to avoid density underflow.

## Source relationships

- The original README identifies Sonjoy Das's FIM calculation as the basis of this work. [estimateFIM.html](../estimateFIM.html) preserves a MATLAB-published source listing with Spall and Das method variants; it is background material, not a required runtime dependency.
- That listing requests citation of Das, Spall, and Ghanem's article identified in [CITATION.cff](../CITATION.cff), [DOI: 10.1016/j.csda.2009.09.018](https://doi.org/10.1016/j.csda.2009.09.018), when using the file. The navigation scripts use its log-likelihood perturbation pattern; they do not expose every method variant in the listing.
- [laprnd.m](../laprnd.m) retains its original Elvis Chen attribution and inverse-transform Laplace sampler.
- The APISAT citation is supported by the original README, script headers, and [author-uploaded paper](https://www.researchgate.net/publication/321825093_Monte-Carlo_Calculation_of_Cramer-Rao_Bound_for_non-Gaussian_Recursive_Filtering).

## Reuse terms

The root [MIT license](../LICENSE) names Youngjoo Kim as copyright holder. The third-party HTML listing and Laplace sampler contain attribution but no separate license grant; their original licensing terms have not been established here, and the root license alone does not establish those upstream permissions.

Retain existing notices when adapting source, and resolve upstream permissions for reuse that relies on those third-party portions. Reference materials can be kept in the ignored `/ref/` directory; none are needed to run the synthetic terrain model.
