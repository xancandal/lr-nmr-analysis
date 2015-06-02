function [y] = pd_Matxxxx(mode,m,n,x,tvec,Tvec)

% [y] = pd_Matxxxx(mode,m,n,x,tvec,Tvec);
%
% Compute the product matrix vector
%
%   y = A * x
%
% where
%   y      is a given m-vector;
%   A      is an m x n matrix;
%   x      is a given n-vector;
%
% such as A is the kernel
%
%   A(s,t) = exp(-t/T)*deltaT                      0 < s,t < 1
%
% of a first Integral Fredholm Kind (IFK) equation discretized by Riemann sums 
%
% INPUT ARGUMENTS:
%   mode             returns  y = A*x (mode=1)  or  y = A'*x (mode=2)
%   m                is the number of data points, i.e. the size of t(i)
%                    or number of echoes. These must be sufficiently close 
%                    each other to capture high decays, and should be enough 
%                    to capture low decays
%   n                discretizacion size of f(T_2), or n - the size of x(j)
%   x                vector (n-size) which is multiplied by the matrix
%   tvec             is the vector of time points
%   Tvec             is the vector of T_2 points
%
% OUTPUT ARGUMENTS:
%   y                is the vector (m-size) result of multiplying the 
%                    matrix A with the vector x

%--------------------------------------------------------------------------
% 13 May 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

dT=diff(Tvec);dT=[dT,dT(end)];         % Delta T_2
if mode==1
    y=zeros(m,1);
    for k=1:m
        a=exp(-tvec(k)./Tvec).*dT;     % matrix unscaled
        y(k)=a*x;
    end
else
    y=zeros(n,1);
    for k=1:n
        a=exp(-tvec./Tvec(k))*dT(k);   % matrix unscaled
        y(k)=a'*x;
    end
end

%--------------------------------------------------------------------------
% End function pd_Matxxxx
%--------------------------------------------------------------------------