function [ L ] = loglikelihood_TRN( Z_est, Z , sig_Z, Dist )
%LOGLIKELIHOOD_TRN 이 함수의 요약 설명 위치
%   자세한 설명 위치

if (Dist == 0) % Gaussian
    var = sig_Z^2;

    p = (1/sqrt(2*pi*var)) * exp( -(Z - Z_est)^2/(2*var) );
else
    b = sig_Z/sqrt(2);
    
    p = (1/(2*b)) * exp( -abs(Z - Z_est)/b);
    
end

L = log(p);
end

