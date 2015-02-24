clear;
%%% Dirac Test (1 point) ..............................................................
t2Values = [1,2];        % T_2 values
amplitudes = [2,1];      % f(T_2) values
snr = 1;                 % Signal-to-noise ratio

nbr_data_points = 1000;  % Discretizacion size of CPMG Signal, i.e. S(t) (m-size)
sampling_period = 0.01;  % Delta t

n = 200;                 % Discretizacion size of amplitudes, i.e. f(T2) (n-size)
deltaT = 0.01;           % Delta T2

% Create Synthetic CPMG Signal
data = CPMGSynthetic_gauss(t2Values, amplitudes, sampling_period, nbr_data_points, n, snr, deltaT);

% Run PDCO
A = zeros(nbr_data_points,n); 
[A,x,z] = PDCOSynthetic_gauss( data, sampling_period, n, snr, deltaT, 0);

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
plot(sampling_period*(0:m-1),A*x,'r');                  % ... PDCO solution
plot(sampling_period*(0:m-1),data(:,2), 'k');           % ... generated solution
hold off;
legend('from exact sol', 'from PDCO sol', 'from Generated sol', 'Location','northeast', 'Orientation','horizontal');
legend('boxoff');
xlabel('t[sec]');
ylabel('Signal S(t)');
