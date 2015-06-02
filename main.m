clc; clear; 

%%% Values of T2 and f(T2) of Signal no. 5 
t2Values = [8.73e-3 21.54e-3 53.15e-3];
amplitudes = [400 400 200];
snr_desired = 18829;

m = 16384;          % size of CPMG Signal, i.e. S(t) size
dt = 1e-4;          % Delta t (or) sampling_period
n = 2^8;            % size of amplitudes, i.e. f(T2) size

%-----------------------------------------------------------------------------

%%% Create the vector of time points and Delta T_2
t = 0:dt:(m-1)*dt;
T =logspace(-4,0,n);

%%% Function handler for operator matrix
A = @(mode,m,n,x) pd_Matxxxx(mode,m,n,x,t,T);

%%% Create Synthetic CPMG Signal without noise and PDFs
[signal,gaussian] = create_signal(t2Values,amplitudes,t,T,n);

%%% Add certain noise level to the Synthetic CPMG Signal for a desired SNR
for noise_level=2e-5:1e-4:2e-1
    [data,snr] = add_noise(signal, t, noise_level);
    if snr < snr_desired
        break;
    end
end
fprintf('\n   Level noise:  %11.4e \n', noise_level)
fprintf('\n   Signal-to-noise ratio:  %11.4e \n', snr)

%%% Run PDCO with scaling
[xx,normr] = pdcoRun( data, n, snr, T);
%% Saving results (Octave)
save xx.mat xx -mat7-binary;
save normr.mat normr -mat7-binary;
%% Saving results (Matlab)
% save xx.mat xx -v7.3;
% save normr.mat normr -v7.3;

%%% Run PDCO with scaling and without regularization L1
[xL2,normrL2] = pdcoRun_L2( data, n, snr, T);
%% Saving results (Octave)
save xL2.mat xL2 -mat7-binary;
save normrL2.mat normrL2 -mat7-binary;
%% Saving results (Matlab)
% save xL2.mat xL2 -v7.3;
% save normrL2.mat normrL2 -v7.3;

%-----------------------------------------------------------------------------

%%% Plot settings
set(0,'defaulttextinterpreter','tex')
set(0,'DefaultTextFontSize', 10)
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultAxesFontSize', 10)
set(0,'DefaultAxesFontName','Times')
set(0, 'defaultlinelinewidth', 2)       % default 1 (Octave), 0.5 (Matlab)

figure
%%% Plot T2 Components and ther Amplitudes building from...
hold off; plot(T,xx,'--b');             % PDCO solution.
hold on;  plot(T,xL2,'-.r');            % L2-PDCO sol.
hold on;  plot(T,sum(gaussian,1),'k');  % generated sol.
axis([1e-4 max(T) 0 450]);
set(gca,'XScale','log');                % plot linear data on logarithmic axes
legend ('CO', 'L2', 'SG',...
  'Location', 'northeast', 'Orientation', 'vertical');
legend ('boxoff');
legend('left');                         % only for Octave
%%% Axis description
xlabel('Relaxation[sec]');
ylabel('Intensity');
%%% To change orientation
orient landscape
%%% Get rid of the white margin in saveas output (for Octave)
dpi = get (0, 'screenpixelsperinch');
pos = get (gcf, 'position');
papersize = pos(3:4)./dpi;
set (gcf, 'papersize', papersize)
set (gcf, 'paperposition', [0, 0, papersize]) 
%%% Saving plot (Octave)
print -landscape -dpdf t2_pdco.pdf
%%% Saving plot (Matlab)
%saveTightFigure('t2_pdco.pdf')

figure
%%% Plot Synthetic CPMG Signal, Components and Gaussian Noise 
hold off; plot(t, data(:,2), 'k');
for k = 1:length(t2Values)
    hold on; plot(t, A(1,m,n,gaussian(k,:)'), '--b');
end
axis([1e-4 max(t) -0.05 14]);
set(gca,'XScale','log');                % plot linear data on logarithmic axes
%%% Axis description
xlabel('Time[sec]');
ylabel('Intensity');
%%% To change orientation
orient landscape
%%% Get rid of the white margin in saveas output (for Octave)
dpi = get (0, 'screenpixelsperinch');
pos = get (gcf, 'position');
papersize = pos(3:4)./dpi;
set (gcf, 'papersize', papersize)
set (gcf, 'paperposition', [0, 0, papersize]) 
%%% Saving plot (Octave)
print -landscape -dpdf components.pdf
%%% Saving plot (Matlab)
%saveTightFigure('components.pdf')

figure
%%% Plot Synthetic CPMG Signal building from ...
hold  off; plot(t,A(1,m,n,xx),'--b');   % PDCO solution.
hold  on ; plot(t,A(1,m,n,xL2),'-.r');  % L2-PDCO sol.
hold  on;  plot(t,data(:,2), 'k');      % generated sol.
axis([1e-4 max(t) -0.05 14]);
set(gca,'XScale','log');                % plot linear data on logarithmic axes
legend ('CO', 'L2', 'SG',...
  'Location', 'northeast', 'Orientation', 'vertical');
legend ('boxoff');
legend('left');                         % only for Octave
%%% Axis description
xlabel('Time[sec]');
ylabel('Intensity');
%%% To change orientation
orient landscape
%%% Get rid of the white margin in saveas output (for Octave)
dpi = get (0, 'screenpixelsperinch');
pos = get (gcf, 'position');
papersize = pos(3:4)./dpi;
set (gcf, 'papersize', papersize)
set (gcf, 'paperposition', [0, 0, papersize]) 
%%% Saving plot (Octave)
print -landscape -dpdf decay.pdf
%%% Saving plot (Matlab)
%saveTightFigure('decay.pdf')

figure
%%% Plot Residuals for each iteration
hold  on;  semilogy(normr,'b');      % PDCO solution.
hold  on;  semilogy(normrL2,'--r');  % L2-PDCO sol.
legend ('CO', 'L2',...
  'Location', 'northeast', 'Orientation', 'vertical');
legend ('boxoff');
legend('left');                          % only for Octave
%%% Axis description
xlabel ('No. Iterations')
ylabel ('||r||')
%%% To change orientation
orient landscape
%%% Get rid of the white margin in saveas output (for Octave)
dpi = get (0, 'screenpixelsperinch');
pos = get (gcf, 'position');
papersize = pos(3:4)./dpi;
set (gcf, 'papersize', papersize)
set (gcf, 'paperposition', [0, 0, papersize]) 
%%% Saving plot (Octave)
print -landscape -dpdf residuals.pdf
%%% Saving plot (Matlab)
%saveTightFigure('residuals.pdf')
