function [y] = pd_Mat(mode,m,n,x,A)

% [y] = pd_Mat(mode,m,n,x,A);
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
%   mode   returns  y = A*x (mode=1)  or  y = A'*x (mode=2)

%--------------------------------------------------------------------------
% 19 March 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

if mode==1
    y=A*x;
else
    y=A'*x;
end

%--------------------------------------------------------------------------
% End function pd_Mat
%--------------------------------------------------------------------------