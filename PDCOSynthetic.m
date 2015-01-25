function PDCOSynthetic( name, scale, wait, lambda_one, lambda_two )

% PDCOSynthetic( 'data_file.txt', 0, 1 );
% It loads data_file.txt and runs it on pdco.m.
%    scale = 0  suppresses scaling (OK if A,b,c are well scaled)
%    scale = 1  requests scaling (default)
%    wait  = 0  prevents pdco from waiting (default)
%    wait  = 1  asks pdco to wait to allow parameters to be reset.
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
%           Ordinary LS rows have d2(i) = 1.
%
% Files such as 'data_file.txt' should be accessible via the
% current path.  Such as,
%    addpath ~/matlab/dataset
% if that directory contains some of the datasets.

%--------------------------------------------------------------------------
% 23 Jan 2015: Example CPMG Synthetic test program for pdco.m,
%              Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

  % Read data from text file, i.e. CPMG (time, amplitudes)
  fid = fopen (name, 'r');
  % load('lp_adlittle');        % data file to test PDCO
  formatSpec = '%f %f';
  sizedata = [2 Inf];
  data = fscanf(fid, formatSpec, sizedata);
  fclose(fid);
  data = data';
  t = data(:,1);                % time m-vector t(i)
  m = length(t);                % Output size
  n = m;                        % Input size
  
  cop   =   @(x) pd_Obj( x );   % operador vector
  b = data (:,2);               % m-vector amplitudes x(i)
  
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
  %         K(s,t) = exp(-s/t)                       0 < s,t < 1
  %
  % Reference: G. M. Wing, "A Primer on Integral Equations of the
  % First Kind", SIAM, 1991; p. 109.
  %
  % Discretized by Quadrature method; in particular, for the midpoint rule.
    
  A = zeros(n,n);             % Initialization
  dt = 1/n;
  ac=0;
  bc=1;
  ds = (bc-ac)*dt;
  ti = dt*((1:n) - 0.5);
  s = ac + ds*((1:n) - 0.5);
  sti = ((1:n)-0.5)*dt;
  
  [T,S] = meshgrid(ti,s);     % Set up matrix
  A = dt*exp(-S./T);
  % for i=1:n
  % A(i,:) = dt*exp(-s(i)*ti.^-1);
  % end
  
  % A should be reasonably well scaled i.e. norm(A,inf)=~1
  fprintf('\n\n   norm(A,inf):  %11.4e', norm(A,inf))
  %------------------------------------------------------------------------

  if nargin < 2                 % scaling (default)
     scale = 1;   
  end
  
  if nargin < 3                 % prevents pdco from waiting (default)
     wait  = 0;
  end
  
  if nargin < 4                 % lambda_one (default)
    lambda_one = 2e-2;          % alfa_one*noiselevel (alfa_one = 0.1, noiselevel = 2e-1)
  end
  
  if nargin < 5                 % lambda_two(default)
    lambda_two = 1e-8;          % gamma = 1e-4
  end
  
  % Upper and lower bounds
  bl      = zeros(n,1);           % x >= 0
  bu      = inf(n,1);
  
  % Scaling ...
  bscale  = norm(b,inf);     bscale  = max(bscale,1);
  %oscale  = norm(bscale*c,inf);     oscale  = max(oscale,1);

  if scale
    b       = b /bscale;
    %c       = (bscale * c) / oscale;
    %fprintf('\n\n  Final b and c scales:  %11.4e     %11.4e', bscale, oscale)
    fprintf('\n\n   Final b scale:  %11.4e', bscale)
  end
  
  % Declare variables as global for the function pd_Obj
  global glambda_one 
  glambda_one = lambda_one;

  % Regularization parameters
  gamma   = sqrt(lambda_two);     % Primal regularization
  delta   = 1;                    % 1e-3 or 1e-4 for LP;  1 for Least squares
  d1      = gamma;                % Can be scalar if D1 = d1*I
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
  options.Method      = 1;        % 1=Chol  2=QR  3=LSMR  4=MINRES
  options.OptTol      = 1e-6;     % 1e-6 is typically small enough;
                                  % 1e-5 may be acceptable also
  options.FeaTol      = 1e-6;     % Typically the same as Featol
  options.wait        = 0;        % 1 Allow options to be reviewed some parameters before solve
  options.MaxIter     = 100;      % Maximum iterations of the primal-dual barrier method
  
  % Parameters for LSMR
  options.LSMRatol1   = 1e-6;     % 1e-3 or 1e-4 sometimes works;
                                  % 1e-8 may be needed for LPs
  options.LSMRatol2   = 1e-6;     % is the smallest value atol is reduced to
  options.LSMRMaxIter = 100;      % Default = 10 (* min(m,n))
  
  % [x,y,z,inform,PDitns,CGitns,time] = ...
  %   pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );  % to test PDCO

  [x,y,z,inform,PDitns,CGitns,time] = ...
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
    disp('Waiting in LPnetlib in case you want to look at the solution')
    keyboard
  end

  % Extract the solution from the output vector x
  j = 0;
  t2Values = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2]; % real values of T_2
  for i=1:length(t2Values)-1
      t2diff(i)=t2Values(i+1)-t2Values(i);
  end
  for i=1:n
      % if x(i) > 1 j = j + 1; xsol(j) = x(i); t2(j) = i * log(A(i,j).^-1); end
      if x(i) > 1
          j = j + 1;
          xsol(j) = x(i)*min(t2diff);
          % tsol = i*(t(2)-t(1));
          t2(j) = j*min(t2diff);
      end
  end
  
  % Plot the first set of data in blue (obtained data set)
  figure
  t2 = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2 0.05];  % real values of T_2
  plot( t2, xsol, 'bo');
  % stem(x ,'bo');
  hold on;
  
  % Plot the second set of data in red (real data set)
  t2Values = [1e-3 2.15e-3 4.64e-3 1e-2 2.15e-2 4.64e-2];
  amplitudes = [1.5 3 2.6 3.3 1.3 5.5];
  plot(t2Values, amplitudes, 'r+');
  
  % Add title and axis labels
  title('T2 Components and their Amplitudes');
  xlabel('T2 Components');
  ylabel('Amplitudes');

  % Add a legend
  legend('PDCO output', 'Synthetic CPMG Signal', 2);
   
  hold off;

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
