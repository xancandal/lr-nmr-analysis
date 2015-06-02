function [signal,gaussian] = ...
    create_signal(t2Values,amplitudes,tvec,Tvec,n,wait)

% [signal,gaussian] = ...
%    create_signal(t2Values,amplitudes,tvec,Tvec,n);
%
% Create Synthetic CPMG Signal without noise.
%
% INPUT ARGUMENTS:
%   t2Values         is a array which contains the max values of T_2
%                    components
%   amplitudes       is a array which contains the max values of Amplitude
%                    (or f(T_2)) components
%   tvec             is the vector of time points
%   Tvec             is the vector of T_2 points
%   n                discretizacion size of f(T_2), or n - the size of x(j)
%   wait = 0         prevents it from waiting (default)
%   wait = 1         asks it to wait to allow parameters to be reset
%
% OUTPUT ARGUMENTS:
%   signal           CPMG signal without noise.
%   gaussian         normal probability density functions

%--------------------------------------------------------------------------
% 13 May 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

if nargin < 6        % prevents pdco from waiting (default)
   wait  = 0;
end

%%% Verification data (only for Matlab)
%testCase = matlab.unittest.TestCase.forInteractiveUse;
%verifyEqual(testCase, size(t2Values), size(amplitudes));

%%% Size of time vector t(i)
m = numel(tvec);

%%% Functions handlers for operator matrix
A = @(mode,m,n,x) pd_Matxxxx(mode,m,n,x,tvec,Tvec); 

%%% Mean values
factor=diff(log10(t2Values));

%%% Standard deviation values
sigma0=3e-3; % standard deviation at peak 21.54e-3
sigma=sigma0*10.^[-factor(1),0,factor(1)];

%%% Initialization
pdfNormal = zeros(length(amplitudes), n);
gaussian = zeros(length(amplitudes), n);

%%% Normal probability density functions
pdfNormal(1,:) = normpdf(Tvec, t2Values(1), sigma(1));
gaussian(1,:) = amplitudes(1) * pdfNormal(1,:)/max(pdfNormal(1,:));
pdfNormal(2,:) = normpdf(Tvec, t2Values(2), sigma(2));
gaussian(2,:) = amplitudes(2) * pdfNormal(2,:)/max(pdfNormal(2,:));
pdfNormal(3,:) = normpdf(Tvec, t2Values(3), sigma(3));
gaussian(3,:) = amplitudes(3) * pdfNormal(3,:)/max(pdfNormal(3,:));
gaussians = sum(gaussian,1);

%%% Compute product between the matrix A and the Laplace transform to 
% obtain the Synthetic CPMG Signal without noise
signal_trasp = zeros(m);
signal_trasp = A(1,m,n,gaussians');
signal = signal_trasp';

if wait
  disp('   Waiting in case you want to look at the solution')
  disp('   To terminate the keyboard mode, type dbcont, and press Enter')
  disp('   To terminate keyboard mode and exit, type dbquit, and press Enter')
  keyboard
end

%--------------------------------------------------------------------------
% End function matrix_signal
%--------------------------------------------------------------------------