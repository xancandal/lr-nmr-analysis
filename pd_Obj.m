function [obj,grad,hess] = pd_Obj(x,lambda_one,n)
    
% Computes the objective value, gradient and diagonal Hessian of the 
% linear function lambda_one * norm(x,1).

%--------------------------------------------------------------------------
% 23 Jan 2015: Xan Candal, University of Santiago de Compostela.
%--------------------------------------------------------------------------

obj  = norm(lambda_one.*x,1);
grad = lambda_one.*sign(x);
hess = zeros(n,1);

%--------------------------------------------------------------------------
% End function pd_Obj
%--------------------------------------------------------------------------