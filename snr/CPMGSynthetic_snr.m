function [data,snr] = CPMGSynthetic_snr( t2Values, amplitudes, sampling_period, nbr_data_points, noise_level, n, deltaT, draw, wait )

% CPMGSynthetic_gauss( t2Values, amplitudes, sampling_period, nbr_data_points, snr );
% It creates a Synthetic CPMG Signal stored into data matrix (time, amplitudes)
%    snr                is the Signal-to-noise ratio per sample, in dB
%    t2Values           is a array which contains the T_2 components
%    amplitudes         is a array which contains the Amplitude (or f(T_2)) components
%    sampling_period    it must be a min distance 1.5 times between the sampling
%                       period and the first value of T_2, i.e. the echo time
%    nbr_data_points    is the number of data points, i.e. the size of t(i) or m  
%                       or number of echoes. These must be sufficiently close each
%                       other to capture high decays, and should be enough to capture low decays
%    noise_level        is the level noise (Gaussian)
%    n                  Discretizacion size of f(T_2), i.e. the size of x(j) or n
%    deltaT             Delta T_2
%    draw = 0           prevents plot figures (default)
%    draw = 1           plot figures
%    wait = 0           prevents CPMGSynthetic from waiting (default)
%    wait = 1           asks CPMGSynthetic to wait to allow parameters to be reset

%--------------------------------------------------------------------------
% 29 Jan 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

  if nargin < 8                 % prevents plot figures (default)
     draw = 1;
  end

  if nargin < 9                 % prevents pdco from waiting (default)
     wait  = 0;
  end
  
%%% Vector of time points
t = 0:sampling_period:(nbr_data_points-1)*sampling_period;

%%% Verification data
% testCase = matlab.unittest.TestCase.forInteractiveUse;
% verifyEqual(testCase, size(t2Values), size(amplitudes));

%%% Create CPMG signal without noise
idy = 1./t2Values;
signal = zeros(1, length(t));
k=0;
for k=1:length(amplitudes)
    signal = signal + ( amplitudes(k) * exp(- idy(k) * t) );
end

%%% Gaussian noise generation
noise = randn(1,length(t)) * noise_level;
signal_plus_noise = signal + noise;

%%% Save CPMG (time, amplitudes) data
data(:,1) = t(:);
%data(:,2) = signal_plus_noise(:);   % With noise
data(:,2) = signal(:);               % Without noise

%--------------------------------------------------------------------------

%%% Moving average of 8 points
x=signal_plus_noise;
windowSize = 8;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,x);

%%% Last 6.25% (1024 points if the total amount is 16384) signal's points
ti = t(length(t)-length(t)*0.0625:end);
yi = y(length(y)-length(y)*0.0625:end);

%%% Line linear regression
[r,m,b] = regression(ti,yi);

%%% Last 6.25% signal's points "minus" Line linear regression
xi = zeros(length(yi), 1);
i = 0;
for i=1:length(ti)
    xi(i) = yi(i)-(m*ti(i)+b);
end

%%% RMS value of the noise
noise_rms = rms(xi);

%%% SNR value
signal_rms= rms(signal_plus_noise);
snr = signal_rms/noise_rms;

%--------------------------------------------------------------------------

%%% Normal probability density function
sigma = 0.1;
Ti = deltaT*(1:n);
pdfNormal = zeros(length(amplitudes), n);
gaussians = zeros(length(amplitudes), n);
i = 0;
for i=1:length(amplitudes)
    pdfNormal(i,:) = normpdf(Ti, t2Values(i), sigma);
    gaussians(i,:) = amplitudes(i) * pdfNormal(i,:)/max(pdfNormal(i,:));
end

if draw
    disp('Plotting figures ...')
    
    %%% Ploting CPMG signal and noise
    plot(t,signal,'b')                                 % Plot signal without noise
    hold on
    plot (t,signal_plus_noise,'r')                     % Plot signal with noise
    grid on
    legend('Original signal','Signal with Noise');
    hold off
    %%% Moving average of 8 points
    figure
    plot(t,x,'r')
    hold on
    plot(t,y,'b')
    grid on
    legend('Input Data','Filtered Data','Location','NorthEast')
    title('Plot of Input and Filtered Data')
    hold off
    %%% Last 6.25% (1024 points if the total amount is 16384) signal's points
    figure
    plot(ti,yi,'b')
    %%% Line linear regression
    figure
    plot(ti, m*ti+b,'g');
    hold on
    plot(ti,yi,'b')
    hold off
    %%% Last 6.25% signal's points "minus" Line linear regression
    figure
    plot(ti,xi,'r')
    
    %%% Plot Synthetic CPMG Signal, Components and Gaussian Noise 
    figure
    subplot(211);
    plot(data(:,1), data(:,2), 'b');
    hold on;
    %plot(data(:,1), noise, 'c');
    for k = 1:length(t2Values)
        plot( data(:,1), (amplitudes(k) * exp(- idy(k) * t) ), 'r');
    end
    hold off;
    title('Synthetic CPMG Signal, Components and Gaussian Noise');
    xlabel('t[sec]');
    ylabel('Echo amplitude S(t)');
    %%% Plot T2 Components and their Amplitudes
    gaussians_sum = sum(gaussians,1);
    subplot(212);
    hold on;
    stem(t2Values, amplitudes,'b');
%     for k = 1:length(t2Values)
%         plot(Ti, gaussians(k,:),'r'); 
%     end
    plot(Ti, gaussians_sum,'r');
    %set(gca,'XScale','log');  % Plot liner data on logarithmic axes
    title('T2 Components and their Amplitudes');
    xlabel('T2[sec]');
    ylabel('Amplitudes');
    hold off;
end

if wait
  disp('Waiting in case you want to look at the solution')
  keyboard
end

