function [x,normr] = pdcoRun(data,n,snr,Tvec,scale,wait)

% [x,normr] = pdcoRun(data,n,snr,Tvec);
%
% Solves optimization problems of the form with scaling matrix A and vector b
% without regularization L1
%
%    minimize    phi(x) + 1/2 norm(D2*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D*r = b,   bl <= x <= bu,   r unconstrained,
%
% where
%    phi(x) is a smooth convex function  defined by function pdObj;
%    A      is an m x n matrix defined by matrix or function pdMat;
%    b      is a given m-vector;
%    D, D2  are positive-definite diagonal matrices defined from d, d2.
%           In particular, d indicates the accuracy required for
%           satisfying each row of Ax = b.
%
% D and D2 (via d and d2) provide primal and dual regularization
% respectively.  They ensure that the primal and dual solutions
% (x,r) and (y,z) are unique and bounded.
%
% A scalar d is equivalent to d = ones(m,1), D = diag(d).
% A scalar d2 is equivalent to d2 = ones(n,1), D2 = diag(d2).
% Set d = 1, d = 1e-4 for least-squares problems with bound constraints.
% The problem is then equivalent to
%
%    minimize    phi(x) + 1/2 norm(d2*x)^2 + 1/2 norm(A*x - b)^2
%    subject to  bl <= x <= bu.
%
% More generally, d and d2 may be m and n vectors containing any positive
% values (preferably not too small, and typically no larger than 1).
% Bigger elements of d and d2 improve the stability of the solver.
%
% INPUT ARGUMENTS:
%   data            read CPMG (time, amplitudes) data.
%   n               size of f(T_2), i.e. the size of x(j) or n
%   snr             Signal-to-noise ratio
%   Tvec            is the vector of T_2 points
%   scale = 0       suppresses scaling (OK if A,b,c are well scaled)
%   scale = 1       requests scaling (default)
%   wait  = 0       prevents pdco from waiting (default)
%   wait  = 1       asks pdco to wait to allow parameters to be reset
%
% OUTPUT ARGUMENTS:
%   x               is the primal solution.
%   normr           is a n array with an estimate of the residual: NORM(b-A*x)
%

%--------------------------------------------------------------------------
% 13 May 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

  if nargin < 5                 % scaling (default)
     scale = 1;   
  end
  
  if nargin < 6                 % prevents pdco from waiting (default)
     wait  = 0;
  end

  %%% Read CPMG (time, amplitudes) data
  t = data(:,1);                % time m-vector t(i)
  m = numel(t);                 % size of CPMG Signal or number of echoes
  b = data (:,2);               % amplitudes without scaling
  
  %%% Create explicit matrix
  A = explicit_matrix(t, Tvec );

  %%% Delta T_2
  dT=diff(Tvec)';dT=[dT;dT(end)];
  ind=find(dT<1e-4);
  dT(ind)=1e-4;

  %%% Scaling matrix A
  % Matrix A should be reasonably well scaled i.e. norm(A,inf)=~1
  Ascale=norm(A,inf);
  fprintf('\n   norm(A,inf):  %11.4e', Ascale)
  if scale
    A = A/Ascale;
    b = b/Ascale;
    disp('   ... scaling using norm(A,inf)')
  end

  %%% Scaling vector b (option 2)
  % Vector b should be reasonably well scaled i.e. norm(b,inf)=~1
  bscale=norm(b,inf);
  fprintf('   norm(b,inf):  %11.4e', bscale)
  if scale
    %A=A;
    b = b/bscale;
    disp('   ... scaling using norm(b,inf)')
  end

  %%% Regularization parameters for LNMR problem
  beta       = norm(b,inf);
  alpha_one  = 10;                             % universal coeficients
  alpha_two  = 5;
  %lambda_one = ((alpha_one*beta)/snr);        % for regularization L1
  lambda_one = 0;                              % without regularization L1
  lambda_two = (alpha_two/snr);                % for regularization L2
  d1         = lambda_one*dT;      
  d2         = sqrt(lambda_two*(dT).^2);       % can be scalar if D2 = d2*I  
  d          = ones(m,1);                      % Least-Squares problem
  fprintf('\n   beta value:  %11.4e', beta)
  fprintf('\n   lambda_one parameter:  %11.4e',   lambda_one)
  fprintf('\n   lambda_two parameter:  %11.4e\n', lambda_two)

  %%% Functions handlers for operator vector and matrix
  cop  =  @(x) pd_Obj(x,d1,n);
  cmat =  @(mode,m,n,x) pd_Mat(mode,m,n,x,A);

  % Upper and lower bounds
  bl = zeros(n,1);                % x >= 0
  bu = inf(n,1);
                 
  if scale
    xsize   = norm(b,inf);        % estimate of norm(x,inf) at solution
    zsize   = 1;                  % estimate of norm(z,inf) at solution
  else
    %xsize   = bscale;
    xsize   = 1;
    zsize   = 1; 
  end

  % Generate an initial guess
  x0      = ones(n,1)*xsize;      % initial x
  y0      = zeros(m,1);           % initial y
  z0      = ones(n,1)*zsize;      % initial z

  options             = pdcoSet;  % option set for the function PDCO
  
  % Grab input options
  options.Print       = 1;        % 1 gives output (default). 0 suppresses it.
  options.Method      = 3;        % 1=Chol  2=QR  3=LSMR  4=MINRES
  options.wait        = 0;        % 1 Allow options to be reviewed some 
                                  % parameters before solve
  options.MaxIter     = 100;      % maximum iterations of the primal-dual 
                                  % barrier method
  
  % Parameters for LSMR
  options.LSMRatol1   = 1e-6;     % 1e-3 or 1e-4 sometimes works;
                                  % 1e-8 may be needed for LPs
  options.LSMRatol2   = 1e-6;     % is the smallest value atol is reduced to
  options.LSMRconlim  = 1e8;      % shuts LSMR down early if its matrix is
                                  % ill-conditioned

  %%% Run PDCO algorithm
  [x,y,z,inform,PDitns,CGitns,time,normr] = ...
    pdco( cop,cmat,b,bl,bu,d2,d,options,x0,y0,z0,xsize,zsize );
  
  %%% Unscaling solution vector x because of scaling vector b (option 2)
  if scale                  
    x  = x*bscale;
  end
  
  if wait
    disp('   Waiting in case you want to look at the solution')
    disp('   To terminate the keyboard mode, type dbcont, and press Enter')
    disp('   To terminate keyboard mode and exit, type dbquit, and press Enter')
    keyboard
  end

  %--------------------------------------------------------------------------
  % End function pdcoRun
  %--------------------------------------------------------------------------