function [f,P1x,t] = fft_analysis_vec(x,frame_rate)
% returns the power sspectrum of a vector
% input: 1D signal (across time)
%        frame rate (int)
% output: frequency vector, power spectrum, time vector
% teresa, 17/03/2023

% Check frequency of movement in that reach
T = 1/frame_rate;
L = numel(x); % Length of signal
t = (0:L-1)*T;           % Time vector
f = frame_rate*(0:(L/2))/L; % frequency vector

% 
Y = fft(x);
P2x = abs(Y/L);
P1x = P2x(1:L/2+1);
P1x(2:end-1) = 2*P1x(2:end-1);

end