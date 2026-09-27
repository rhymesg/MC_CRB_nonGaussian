function [ h ] = meas_TRN_1d( x )
%MEAS_TRN_1D Synthetic sinusoidal terrain height for position x.
% Model and units: docs/algorithm.md#measurement-and-likelihood-contracts

    Off = 70;
    Amp = 12;
    Len = 80;

    
    h = Off + Amp*sin(x/Len * 2*pi);
end

