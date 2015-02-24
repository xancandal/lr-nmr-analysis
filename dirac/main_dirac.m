clear;
%%% Dirac Test (1 point) ..............................................................
t2Values = [1]; 
amplitudes = [2];
nbr_data_points = 100;  % Discretizacion size of CPMG Signal (m-size)
sampling_period = 0.1;  % Delta t
deltaT = 0.1;           % Delta T2
n = 20;                 % Discretizacion size of f(T2) (n-size)
snr = 1;                % Signal-to-noise ratio

% Create Synthetic CPMG Signal
data = CPMGSynthetic_dirac(t2Values, amplitudes, sampling_period, nbr_data_points, snr);

% Run PDCO
A = zeros(nbr_data_points,n); 
[A,x,z] = PDCOSynthetic_dirac( data, sampling_period, n, snr, 0);

figure
%%% Plot T2 Components and their Amplitudes
subplot(411);
stem(t2Values, amplitudes, 'r+');
xlabel('T2[sec]');
ylabel('Generated f(T2)');

%%% Plot T2 Components and computed solutions with PDCO
subplot(412);
plot(sampling_period*(1:n),x);
xlabel('T2[sec]');
ylabel('PDCO f(T2)');

%%% Plot T2 Components and exact solutions
subplot(413);
x_exacto = zeros(n,1);
Ti= deltaT*(1:n);
x_exacto(Ti==1)=2;
plot(sampling_period*(1:n),x_exacto/deltaT);
norm(A*x_exacto-data(:,2),inf)
xlabel('T2[sec]');
ylabel('Exact f(T2)');

%%% Plot Synthetic CPMG Signal building from ...
subplot(414);
m= numel(data(:,2));
hold on;
plot(sampling_period*(0:m-1),A*x_exacto/deltaT, 'b');   % ... exact solution
plot(sampling_period*(0:m-1),A*x,'r');                  % ... PDCO solution
plot(sampling_period*(0:m-1),data(:,2), 'k');           % ... generated solution
hold off;
legend('from exact sol', 'from PDCO sol', 'from Generated sol', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
xlabel('t[sec]');
ylabel('Signal S(t)');
