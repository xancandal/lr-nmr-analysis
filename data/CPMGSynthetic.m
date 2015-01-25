%%% it must be a distance 1.5 times between sampling rate and the first value of t2
sampling_rate = 1e-5;
%%% number of data points, these must be sufficiently close each other to capture
%%% high decays, and should be enough to capture low decays
nbr_data_points = 2^10;
%%% level noise
noise_level = 2e-1;
%%% vector of time points
t = 0:sampling_rate:(nbr_data_points-1)*sampling_rate;
t2Values = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2];
% t2Values = [.8e-3 2.20e-3 2.55e-3 1e-2 2.15e-2 4.64e-2];
amplitudes = [1.5 3 2.6 3.3 1.3 5.5];
% amplitudes = [2 3 4 2 3 4];
testCase = matlab.unittest.TestCase.forInteractiveUse;
%%% Verification
verifyEqual(testCase, size(t2Values), size(amplitudes));
idy = 1./t2Values;
signal = zeros(1, length(t));
for k=1:length(amplitudes)
    signal = signal + ( amplitudes(k) * exp(- idy(k) * t) );
end
%%% Gaussian noise generation
noise = randn(1,length(t)) * noise_level;
data(:,1) = t(:);
signal2 = signal + noise;
data(:,2) = signal2(:);
clf('reset');

subplot(211);
plot(data(:,1), data(:,2), 'b');
hold on;
plot(data(:,1), noise, 'c');
for k = 1:length(t2Values)
    plot( data(:,1), (amplitudes(k) * exp(- idy(k) * t) ), 'r');
end
hold off;
title('Synthetic CPMG Signal, Components and Gaussian Noise');

subplot(212);
semilogx(t2Values, amplitudes, 'or','markersize',10);
axis([0.5*min(t2Values) 2*max(t2Values) 0.8*min(amplitudes) 1.2*max(amplitudes)]);
title('T2 Components and their Amplitudes');

% Write data to text file (time, amplitudes)
data_output = data';
fid = fopen ('data_file.txt', 'w');
fprintf(fid,'%.15f %.15f\n', data_output);
fclose(fid);
