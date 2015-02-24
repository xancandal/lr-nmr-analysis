% Dados: noise_level
% t2Values, amplitudes
% sampling_period, nbr_data_points
%
% Devuelve: SNR
% -------------------------------------------------------------------------

t2Values = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2];
amplitudes = [1.5 3 2.6 3.3 1.3 5.5];

sampling_period = 1e-5;
nbr_data_points = 2^10;

t = 0:sampling_period:(nbr_data_points-1)*sampling_period;

idy = 1./t2Values;
signal = zeros(1, length(t));
k=0;
for k=1:length(amplitudes)
    signal = signal + ( amplitudes(k) * exp(- idy(k) * t) );
end

%%% Level noise
noise_level = 2e-1;
%%% Gaussian noise generation
noise = randn(1,length(t)) * noise_level;
signal_plus_noise = signal + noise;  % With noise
%data(:,2) = signal(:);              % Without noise

%%% Ploting CPMG signal and noise
plot(t,signal,'b')                                 % Plot signal without noise
hold on
plot (t,signal_plus_noise,'r')                     % Plot signal with noise
grid on
legend('Original signal','Signal with Noise');
hold off

%%% Save CPMG (time, amplitudes) data
data(:,1) = t(:);
data(:,2) = signal_plus_noise(:);

%%% Moving average of 8 points
x=signal_plus_noise;
windowSize = 8;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,x);

figure
plot(t,x,'r')
hold on
plot(t,y,'b')

grid on
legend('Input Data','Filtered Data','Location','NorthEast')
title('Plot of Input and Filtered Data')
hold off

%%% Last 6.25% (1024 points if the total amount is 16384) signal's points
ti = t(length(t)-length(t)*0.0625:end);
yi = y(length(y)-length(y)*0.0625:end);
figure
plot(ti,yi,'b')

%%% Line linear regression
[r,m,b] = regression(ti,yi);

figure
plot(ti, m*ti+b,'g');
hold on
plot(ti,yi,'b')
hold off

%%% Last 6.25% signal's points "minus" Line linear regression
xi = zeros(length(yi), 1);
i = 0;
for i=1:length(ti)
    xi(i) = yi(i)-(m*ti(i)+b);
end

figure
plot(ti,xi,'r')

%%% RMS value of the noise
noise_rms = rms(xi);

%%% SNR value
signal_rms= rms(signal_plus_noise);
snr = signal_rms/noise_rms
