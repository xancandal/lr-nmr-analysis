clear; 
%%% Test2 (3 points) ..................................................................
t2Values = [1,0.25,2];     % T_2 values
amplitudes = [2,0.5,1];      % f(T_2) values

snr_desired = 1200;             % Desired Signal-to-noise

nbr_data_points = 1000;   % Discretizacion size of CPMG Signal, i.e. S(t) (m-size)
sampling_period = 0.01;   % Delta t
n = 200;                  % Discretizacion size of amplitudes, i.e. f(T2) (n-size)
deltaT = 0.01;            % Delta T2

%%% Create Synthetic CPMG Signal for a desired SNR
for noise_level=2e-5:1e-4:2e-1
    [data,snr] = CPMGSynthetic_test(t2Values, amplitudes, sampling_period, nbr_data_points, noise_level, n, deltaT, 0);
    if snr < snr_desired
        break;
    end
end
fprintf('\n   Level noise:  %11.4e \n', noise_level)
fprintf('\n   Signal-to-noise ratio:  %11.4e \n', snr)

% Run PDCO with scaling
A = zeros(nbr_data_points,n); 
[A,xs] = PDCOSynthetic_test( data, sampling_period, n, snr, deltaT);
% Run PDCO without scaling
B = zeros(nbr_data_points,n); 
[B,x] = PDCOSynthetic_test( data, sampling_period, n, snr, deltaT, 0);

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
gaussians_sum = sum(gaussians,1);

figure
%%% Plot T2 Components and their Amplitudes
subplot(411);
hold on;
stem(t2Values, amplitudes,'b');
% for k = 1:length(t2Values)
%    plot(Ti, gaussians(k,:),'r'); 
% end
plot(Ti, gaussians_sum,'r');
%set(gca,'XScale','log');  % Plot liner data on logarithmic axes
title('T2 Components and their Amplitudes');
xlabel('T2[sec]');
ylabel('Generated f(T2)');
hold off;

%%% Plot T2 Components and computed solutions with PDCO
subplot(412);
hold on;
plot(sampling_period*(1:n),xs,'r');
plot(sampling_period*(1:n),x,'b');
hold off;
legend('with scaling', 'without scaling', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
xlabel('T2[sec]');
ylabel('PDCO f(T2)');

%%% Plot T2 Components and exact solutions
subplot(413);
x_exacto = zeros(n,length(t2Values));
Ti = deltaT*(1:n);
k = 0;
hold on;
for k=1:length(t2Values)
    x_exacto(Ti==t2Values(k),k)=amplitudes(k);
    plot(sampling_period*(1:n),x_exacto(:,k)/deltaT);
    fprintf('\n   Error of x(%d):  %11.4e', k, norm(A*x_exacto(:,k)-data(:,2),inf) );
end
fprintf('\n');
xlabel('T2[sec]');
ylabel('Exact f(T2)');
hold off;

%%% Plot Synthetic CPMG Signal building from ...
subplot(414);
m = numel(data(:,2));
x_exacto_sum = sum(x_exacto,2); 
hold on;
plot(sampling_period*(0:m-1),A*x_exacto_sum/deltaT, 'b');   % ... exact solution
plot(sampling_period*(0:m-1),A*x,'r');                      % ... PDCO solution
plot(sampling_period*(0:m-1),A*xs,'g');                     % ... PDCO scaled solution
plot(sampling_period*(0:m-1),data(:,2), 'k');               % ... generated solution
hold off;
set(gca,'XScale','log');  % Plot liner data on logarithmic axes
legend('from exact sol', 'from PDCO sol', 'from PDCO scaled sol', 'from Generated sol', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
xlabel('t[sec]');
ylabel('Signal S(t)');
