function [A,x,normA_out] = PDCOSynthetic_signal3( data, sampling_period, n, snr, deltaT, scale, wait)

% PDCOSynthetic( data, sampling_period, n, snr );
% It loads data matrix and runs it on pdco.m.
%    sampling_period    it must be a min distance 1.5 times between the sampling 
%                       period and the first value of T_2, i.e. the echo time
%    n                  Discretizacion size of f(T_2), i.e. the size of x(j) or n
%    SNR                Signal-to-noise ratio
%    deltaT             Delta T_2
%    scale = 0          suppresses scaling (OK if A,b,c are well scaled)
%    scale = 1          requests scaling (default)
%    wait  = 0          prevents pdco from waiting (default)
%    wait  = 1          asks pdco to wait to allow parameters to be reset
%
% PDCO solves optimization problems of the form
%
%    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D2*r = b,   0 <= x <= inf,   r unconstrained,
%
% where
%    phi(x) is a smooth convex function  defined by function pdObj;
%    A      is an m x n matrix defined by matrix;
%    b      is a given m-vector loaded from data_file.txt;
%    D1, D2 are positive-definite diagonal matrices defined from d1, d2.
%

%--------------------------------------------------------------------------
% 23 Jan 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

  % Read CPMG (time, amplitudes) data
  % load('lp_adlittle');        % data file to test PDCO
  t = data(:,1);                % time m-vector t(i)
  m = length(t);                

  cop   =   @(x) pd_Obj( x );   % Operador vector
  b = data (:,2);               
  
  % A       = Problem.A;        % datas to test PDCO
  % b       = Problem.b;
  % bl      = Problem.aux.lo;
  % bu      = Problem.aux.hi;   
  % c       = Problem.aux.c;
  % [m,n]   = size(A);
  
  %------------------------------------------------------------------------
  % Discretization of a first kind Fredholm integral equation with
  % kernel K given by
  %
  %         K(s,t) = exp(-t/T)                       0 < s,t < 1
  %
  % Discretized by Riemann sums 
  % (or)
  % Discretized by Quadrature method; in particular, for the midpoint rule. 
  % Reference: G. M. Wing, "A Primer on Integral Equations of the
  % First Kind", SIAM, 1991. 
  %------------------------------------------------------------------------
 
  C = zeros(m,n);                      % Initialization              
  dT = deltaT;                         % Delta T 
  dt = sampling_period;                % Delta t
  Ti = dT*(1:n);
  ti = (0:m-1)*dt;
  % dT = sampling_period;              % Midpoint rule (alternative)
  % ac=0;
  % bc=1;
  % dt = (bc-ac)*dT;
  % 
  % Ti = dT*((1:n) - 0.5);
  % ti = ac + dt*((0:m-1) - 0.5);
  
  [Tmat,tmat] = meshgrid(Ti,ti);     % Set up matrix
  C = dT*exp(-tmat./Tmat);
  
  % A should be reasonably well scaled i.e. norm(A,inf)=~1
  fprintf('\n   norm(A,inf):  %11.4e', norm(C,inf))
  % Condition number of A
  %cond(C);
  
  % Remove the elements of A < eps(precision)
  Ad = zeros(m,n); 
  for i=1:m
    for j=1:n
        if C(i,j)>eps('double')
        %if C(i,j)>eps('single')
            Ad(i,j) = C(i,j);
        end
    end
  end
  
  % Converts a full matrix to sparse form by squeezing out any zero elements  
  A = sparse(Ad);
  
  %------------------------------------------------------------------------
  
  if nargin < 6                 % scaling (default)
     scale = 1;   
  end
  
  if nargin < 7                 % prevents pdco from waiting (default)
     wait  = 0;
  end
  
  % Upper and lower bounds
  bl      = zeros(n,1);           % x >= 0
  bu      = inf(n,1);
  
  % Scaling ...
  bscale  = norm(b,inf);     bscale  = max(bscale,1);
  fprintf('\n   norm(b,inf):  %11.4e', bscale)
  %oscale  = norm(bscale*c,inf);     oscale  = max(oscale,1);
  %fprintf('\n\n  Final b and c scales:  %11.4e     %11.4e\n', bscale, oscale)

  if scale
    b       = b /bscale;
    disp('   ... scaling vector b')
    %c       = (bscale * c) / oscale;
  end
  
  % Regularization parameters
  beta       = norm(b,inf);
  alpha_one  = 10;            % universal coeficients
  alpha_two  = 5;
  lambda_one = (alpha_one*beta)/snr;
  lambda_two = alpha_two/snr;
  %lambda_one = 0;            % Without regularization L1
  %lambda_two = 0;            % Without regularization L2
  fprintf('\n   beta value:  %11.4e', beta)
  fprintf('\n   lambda_one parameter:  %11.4e', lambda_one)
  fprintf('\n   lambda_two parameter:  %11.4e', lambda_two)
  
  % Declare variables as global for the function pd_Obj
  global glambda_one 
  glambda_one = lambda_one;

  % Regularization parameters
  % Typically, d1 = d2 = 1e-4 (preferably not too small, and typically no larger than 1).
  % Set d1 = 1e-4, d2 = 1 for least-squares problems with bound constraints.
  gamma   = sqrt(lambda_two); 
  if gamma > 1e-4
      gamma = 1e-4;               % lambda_two = 1e-8
  end
  d1      = gamma;                % Can be scalar if D1 = d1*I
  fprintf('\n   Final d1 PDCO parameter:  %11.4e\n', d1)
  delta   = 1;                    % Least-Squares problem
  d2      = delta*ones(m,1);      
  
  % Scaling ...
  if scale
      d1 = d1*bscale;
      %d1 = d1/sqrt(oscale);
      d2 = d2/bscale;
      %d2 = d2*sqrt(oscale);
  end
      
  % Generate an initial guess
  x0      = zeros(n,1);           % Initial x
  x0      = max(x0,bl);
  x0      = min(x0,bu);
  y0      = zeros(m,1);           % Initial y
  z0      = zeros(n,1);           % Initial z
  
  if scale
    xsize   = 1;               % Estimate of norm(x,inf) at solution
    zsize   = 1;               % Estimate of norm(z,inf) at solution
  else
    xsize   = bscale;
    zsize   = 1; 
    %zsize   = oscale;
  end

  options             = pdcoSet;  % Option set for the function PDCO
  
  % Grab input options
  options.mu0         = 1e-0;     % Initial mu (ABSOLUTE VALUE) for solving scaled problem
  options.Method      = 3;        % 1=Chol  2=QR  3=LSMR  4=MINRES
  options.OptTol      = 1e-6;     % 1e-6 is typically small enough;
                                  % 1e-5 may be acceptable also
  options.FeaTol      = 1e-6;     % Typically the same as Featol
  options.wait        = 0;        % 1 Allow options to be reviewed some parameters before solve
  options.MaxIter     = 50;      % Maximum iterations of the primal-dual barrier method
  
  % Parameters for LSMR
  options.LSMRatol1   = 1e-6;     % 1e-3 or 1e-4 sometimes works;
                                  % 1e-8 may be needed for LPs
  options.LSMRatol2   = 1e-6;     % is the smallest value atol is reduced to
  options.LSMRMaxIter = 100;      % Default = 10 (* min(m,n))
  
  % [x,y,z,inform,PDitns,CGitns,time] = ...
  %   pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );  % to test PDCO

  [x,y,z,inform,PDitns,CGitns,time,normA_out] = ...
    pdco( cop,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );
  
  if scale                  % Unscale b,x,y,z
    b  = b *bscale;
    %c  = c *oscale;
    x  = x *bscale;
    y  = y /bscale;
    %y  = y *oscale;
    z  = z /bscale;
    %z  = z *oscale;
  end
  
  if wait
    disp('   Waiting in case you want to look at the solution')
    keyboard
  end

%--------------------------------------------------------------------------
% End function pdcotest
%--------------------------------------------------------------------------

function [obj,grad,hess] = pd_Obj( x )
%        
% computes the objective value, gradient and diagonal Hessian
% of the linear function lambda_one * norm(x,1), where lambda is stored in
% global glambda_one.
% For use with pdcotest.m (which uses pdsco.m).

%--------------------------------------------------------------------------
% 23 Jan 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

global glambda_one

lambda_one = glambda_one;

n    = length(x);
obj  = lambda_one * norm(x,1);
% obj  = lambda_one * sum(x);
grad = lambda_one * ones(n,1);
hess = zeros(n,1);

%--------------------------------------------------------------------------
% End function pd_Obj
%--------------------------------------------------------------------------
