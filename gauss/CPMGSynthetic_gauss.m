function [data] = CPMGSynthetic_gauss( t2Values, amplitudes, sampling_period, nbr_data_points, n, snr, deltaT, draw, wait )

% CPMGSynthetic_gauss( t2Values, amplitudes, sampling_period, nbr_data_points, snr );
% It creates a Synthetic CPMG Signal.
%    t2Values           is a array which contains the T_2 components
%    amplitudes         is a array which contains the Amplitude (or f(T_2)) components
%    sampling_period    it must be a min distance 1.5 times between the sampling
%                       period and the first value of T_2, i.e. the echo time
%    nbr_data_points    is the number of data points, i.e. the size of t(i) or m  
%                       or number of echoes. These must be sufficiently close each
%                       other to capture high decays, and should be enough to capture low decays
%    n                  Discretizacion size of f(T_2), i.e. the size of x(j) or n
%    SNR                is the Signal-to-noise ratio per sample, in dB
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
testCase = matlab.unittest.TestCase.forInteractiveUse;
verifyEqual(testCase, size(t2Values), size(amplitudes));

%%% Create CPMG signal without noise
idy = 1./t2Values;
signal = zeros(1, length(t));
k=0;
for k=1:length(amplitudes)
    signal = signal + ( amplitudes(k) * exp(- idy(k) * t) );
end

%%% Gaussian noise generation
% noise = randn(1,length(t)) * noise_level;
% signal_plus_noise = signal + noise;
% clf('reset');

%%% Add white Gaussian noise with SNR, which is measured as a ratio, asuming signal with power 0 dBW
signal_plus_noise = awgn(signal,snr,'measured','linear');

%%% Save CPMG (time, amplitudes) data
data(:,1) = t(:);
%data(:,2) = signal_plus_noise(:);  % With noise
data(:,2) = signal(:);              % Without noise

% plot(t,signal,t,signal_plus_noise)             % Plot both signals.
% legend('Original signal','Signal with AWGN');
% plot (data(:,1),data(:,2))                     % Plot Orignal signal
% legend('Original signal');

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
    clf('reset');
    
    %%% Plot Synthetic CPMG Signal, Components and Gaussian Noise 
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

