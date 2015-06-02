function [data,snr] = add_noise(signal,t,noise_level,draw,wait)

% [data,snr] = add_noise(signal,t,noise_level);
%
% It creates a Synthetic CPMG Signal stored into a matrix (time,amplitudes)
%
% INPUT ARGUMENTS:
%    signal             CPMG signal without noise.
%    t                  is the vector of time points
%    noise_level        is the level noise (Gaussian)
%    draw = 0           prevents plot figures (default)
%    draw = 1           plot figures
%    wait = 0           prevents it from waiting (default)
%    wait = 1           asks it to wait to allow parameters to be reset
%
% OUTPUT ARGUMENTS:
%    data               it loads data matrix.
%    snr                is the Signal-to-noise ratio per sample, in dB

%--------------------------------------------------------------------------
% 5 March 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

  if nargin < 4         % prevents plot figures (default)
     draw = 0;
  end

  if nargin < 5         % prevents pdco from waiting (default)
     wait  = 0;
  end
  
%%% Gaussian noise generation
noise = randn(1,length(t)) * noise_level;
signal_plus_noise = signal + noise;

%%% Save CPMG (time, amplitudes) data
data(:,1) = t(:);
data(:,2) = signal_plus_noise(:);     % with noise
%data(:,2) = signal(:);               % without noise

%--------------------------------------------------------------------------

%%% Moving average of 8 points
x=signal_plus_noise;
windowSize = 8;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,x);

%%% Last 1024 signal's points (6.25% if the total amount is 16384)
% temp=(nbr_data_points/1024)^-1;
% ti = t(round(length(t)-length(t)*temp):end);
% yi = y(round(length(y)-length(y)*temp):end);

%%% Last 6.25% (1024 points if the total amount is 16384) signal's points
temp = 0.0625;
ti = t(round(length(t)-length(t)*temp):end);
yi = y(round(length(y)-length(y)*temp):end);

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

if draw
    disp('Plotting figures ...')
    
    %%% Ploting CPMG signal and noise
    plot(t,signal,'b')                  % plot signal without noise
    hold on
    plot (t,signal_plus_noise,'r')      % plot signal with noise
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
    %%% Last 6.25%  signal's points (1024 points if total amount is 16384)
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
end

if wait
  disp('   Waiting in case you want to look at the solution')
  disp('   To terminate the keyboard mode, type dbcont, and press Enter')
  disp('   To terminate keyboard mode and exit, type dbquit, and press Enter')
  keyboard
end

%--------------------------------------------------------------------------
% End function add_noise
%--------------------------------------------------------------------------
