   function [A] = explicit_matrix(tvec,Tvec)
  %
  % Discretization of a first Integral Fredholm Kind (IFK) equation with
  % kernel K given by
  %
  % A(s,t) = exp(-t/T)*dT                         0 < s,t < 1
  %
  % Discretized by Riemann sums
  
  %--------------------------------------------------------------------------
  % 13 May 2015: Xan Candal, University of Santiago de Compostela.
  %--------------------------------------------------------------------------

  dT=diff(Tvec);dT=[dT,dT(end)];      % Delta T_2

  [Tmat,tmat] = meshgrid(Tvec,tvec);     
  [dTmat,tmat] = meshgrid(dT,tvec);    
  A = dTmat.*exp(-tmat./Tmat);        % matrix unscaled

  %--------------------------------------------------------------------------
  % End function explicit_matrix
  %--------------------------------------------------------------------------