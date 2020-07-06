function [ h ] = meas_TRN_1d( x )
%MEAS_TRN_1D 이 함수의 요약 설명 위치
%   sine function

    Off = 70;
    Amp = 12;
    Len = 80;

    
    h = Off + Amp*sin(x/Len * 2*pi);
end

