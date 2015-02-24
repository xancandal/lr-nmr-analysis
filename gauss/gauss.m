%t2Values = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2];
%amplitudes = [1.5 3 2.6 3.3 1.3 5.5];
t2Values = [1,0.5,0.25,2];        % T_2 values
amplitudes = [2,1,0.5,1];      % f(T_2) values

noise_level = 2e-1;       % Level noise

% nbr_data_points = 2^10;
% sampling_period = 1e-5
nbr_data_points = 1000;  % Discretizacion size of CPMG Signal, i.e. S(t) (m-size)
sampling_period = 0.01;  % Delta t

n = 200;                  % Discretizacion size of amplitudes, i.e. f(T2) (n-size)
deltaT = 0.01;            % Delta T2

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

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
% for k = 1:length(t2Values)
%    plot(Ti, gaussians(k,:),'r'); 
% end
plot(Ti, gaussians_sum,'r');
%set(gca,'XScale','log');  % Plot liner data on logarithmic axes
title('T2 Components and their Amplitudes');
xlabel('T2[sec]');
ylabel('Amplitudes');
hold off;