function [ L ] = loglikelihood_TRN( Z_est, Z , sig_Z, Dist )
%LOGLIKELIHOOD_TRN Scalar Gaussian or Laplace observation log density.
% Inputs and numerical limitations: docs/algorithm.md#measurement-and-likelihood-contracts

if (Dist == 0) % Gaussian
    var = sig_Z^2;

    p = (1/sqrt(2*pi*var)) * exp( -(Z - Z_est)^2/(2*var) );
else
    b = sig_Z/sqrt(2);
    
    p = (1/(2*b)) * exp( -abs(Z - Z_est)/b);
    
end

L = log(p);
end

