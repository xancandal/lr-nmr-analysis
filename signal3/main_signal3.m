clc; clear; 
% Working with single precision instead of default double precision
%feature('SetPrecision', 24);
%%% Signal 3 ..............................................................
t2Values = [3.54e-3 21.54e-3 131.11e-3];
amplitudes = [400 400 200];
snr_desired = 17395;
% snr_2 = 1722;
% snr_3 = 583;
% snr_4 = 285;

%%% it must be a min distance 1.5 times between sampling rate and the first value of t2
sampling_rate = 10^(log10(min(t2Values))-1.5);  % 1.5 times

nbr_data_points = 16384;         % Discretizacion size of CPMG Signal, i.e. S(t) (m-size)
sampling_period = 1e-4;         %% Delta t
n = 2560;                       %% Discretizacion size of amplitudes, i.e. f(T2) (n-size)
deltaT = sampling_period;        % Delta T2

%%% Create Synthetic CPMG Signal for a desired SNR
for noise_level=2e-5:1e-4:2e-1
    [data,snr] = CPMGSynthetic_signal3(t2Values, amplitudes, sampling_period, nbr_data_points, noise_level, n, deltaT, 0);
    if snr < snr_desired
        break;
    end
end
%[data,snr] = CPMGSynthetic_signal3(t2Values, amplitudes, sampling_period, nbr_data_points, 0.2, n, deltaT, 1); % SNR 200
fprintf('\n   Level noise:  %11.4e \n', noise_level)
fprintf('\n   Signal-to-noise ratio:  %11.4e \n', snr)

% Run PDCO with scaling
As = zeros(nbr_data_points,n); 
[As,xs,normA] = PDCOSynthetic_signal3( data, sampling_period, n, snr, deltaT);
% Run PDCO without scaling
AA = zeros(nbr_data_points,n);
[AA,x] = PDCOSynthetic_signal3( data, sampling_period, n, snr, deltaT, 0);

%--------------------------------------------------------------------------

%%% Normal probability density function
%sigma = 0.5;
Ti = deltaT*(1:n);
pdfNormal = zeros(length(amplitudes), n);
gaussians = zeros(length(amplitudes), n);
% i = 0;
% for i=1:length(amplitudes)
%     pdfNormal(i,:) = normpdf(Ti, t2Values(i), sigma);
%     gaussians(i,:) = amplitudes(i) * pdfNormal(i,:)/max(pdfNormal(i,:));
% end
pdfNormal(1,:) = normpdf(Ti, t2Values(1), 0.001);
gaussians(1,:) = amplitudes(1) * pdfNormal(1,:)/max(pdfNormal(1,:));
pdfNormal(2,:) = normpdf(Ti, t2Values(2), 0.005);
gaussians(2,:) = amplitudes(2) * pdfNormal(2,:)/max(pdfNormal(2,:));
pdfNormal(3,:) = normpdf(Ti, t2Values(3), 0.05);
gaussians(3,:) = amplitudes(3) * pdfNormal(3,:)/max(pdfNormal(3,:));
%gaussians_sum = sum(gaussians,1);

% Eigenvalues of A^T*A (i.e. the squares of the singular values of A)
singlevalues = svd(full(As));

%--------------------------------------------------------------------------

figure
%%% Plot T2 Components and their Amplitudes
subplot(211);
hold on;
stem(t2Values, amplitudes,'b');
for k = 1:length(t2Values)
   plot(Ti, gaussians(k,:),'r'); 
end
%plot(Ti, gaussians_sum,'r');
stem(sampling_rate, max(amplitudes), 'g');            % Plot sampling_rate teorical
stem(sampling_period, max(amplitudes), 'y');          % Plot sampling_rate real
axis([0.5*min(sampling_period) 2*max(t2Values) 0.8*min(amplitudes) 1.2*max(amplitudes)]);
set(gca,'XScale','log');                              % Plot liner data on logarithmic axes
title('T2 Components and their Amplitudes');
xlabel('T2[sec]');
ylabel('Generated f(T2)');
hold off;

%%% Plot T2 Components and computed solutions with PDCO
subplot(212);
hold on;
plot(sampling_period*(1:n),xs,'r');                     % ... PDCO scaled solution
plot(sampling_period*(1:n),x,'b');                      % ... PDCO solution
stem(sampling_rate, max(amplitudes), 'g');              % Plot sampling_rate teorical
stem(sampling_period, max(amplitudes), 'y');            % Plot sampling_rate real
axis([0.5*min(sampling_period) 2*max(t2Values) 0.8*min(xs) 1.2*max(xs)]);
set(gca,'XScale','log');                                % Plot liner data on logarithmic axes
legend('with scaling', 'without scaling', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
hold off;
xlabel('T2[sec]');
ylabel('PDCO f(T2)');
% print -djpg t2_pdco.jpg  % Print plot to jpg (Octave)

idy = 1./t2Values;
t = 0:sampling_period:(nbr_data_points-1)*sampling_period;
figure
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

%%% Plot Synthetic CPMG Signal building from ...
subplot(212);
m = numel(data(:,2));
hold on;
plot(sampling_period*(0:m-1),AA*x,'r');                      % ... PDCO solution
plot(sampling_period*(0:m-1),As*xs,'g');                     % ... PDCO scaled solution
plot(sampling_period*(0:m-1),data(:,2), 'k');                % ... generated solution
set(gca,'XScale','log');                                     % Plot linear data on logarithmic axes
legend('from PDCO sol', 'from PDCO scaled sol', 'from Generated sol', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
hold off;
xlabel('t[sec]');
ylabel('Signal S(t)');
% print -djpg decay.jpg  % Print plot to jpg (Octave)

figure
%%% Plot Residuals for each iteration
subplot(211);
plot(log10(normA)  ,'b-')
xlabel('No. of iterations (k)')
ylabel('Residuals norm(Ar_{k-1})')
title('Residuals for each iteration')

%%% Plot the eigenvalues of A^T*A (i.e. the squares of the singular values of A)
subplot(212);
plot(singlevalues.^2,'b.')
xlabel('Eigenvalue number')
ylabel('\lambda(A^T*A)')
title('Eigenvalues of A^T*A')
% print -djpg res_eig.jpg  % Print plot to jpg (Octave)

%%% Plots the sparsity pattern of any matrix A
figure
spy(As);
% print -djpg spy.jpg  % Print plot to jpg (Octave)

%--------------------------------------------------------------------------

%%% Saving results (Octave)
% save As.mat As -mat7-binary;
% save AA.mat AA -mat7-binary;
% save xs.mat xs -mat7-binary;
% save x.mat x -mat7-binary;
