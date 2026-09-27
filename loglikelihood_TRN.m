function [ L ] = loglikelihood_TRN( Z_est, Z , sig_Z, Dist )
%LOGLIKELIHOOD_TRN Scalar Gaussian or Laplace observation log density.
% Inputs and numerical limitations: docs/algorithm.md#measurement-and-likelihood-contracts

if (Dist == 0) % Gaussian
    var = sig_Z^2;

    L = -0.5*log(2*pi*var) - (Z - Z_est)^2/(2*var);
else
    b = sig_Z/sqrt(2);
    
    L = -log(2*b) - abs(Z - Z_est)/b;
    
end

end

