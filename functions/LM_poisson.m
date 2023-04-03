function [a, iter] = LM_poisson(f, Jac, x, v0, tol, maxiter, lambda, factor)
    a = v0'; %IMPORTANT: a must be a column vector!
    n = length(a);
    iter = 0;

    %varargs:
    if nargin < 8
        factor = 10; 
    end
    if nargin < 7
        lambda = 1e-2;
    end
    if nargin < 6                    
        maxiter = 100;
    end
    if nargin < 5 
        tol = 1e-7;
    end
    
    while true
        f_a = f(a);
        
        J = Jacobian_Alex(f,a); 
        %[J, ~,~] = Adigator_f(a);
        %J = Jac(f_a, a); 

        tmp = x(:) ./ f_a(:).^2;
        alpha = transpose(J) * (J .* repmat(tmp,1,5)); %
        alpha = alpha .* (1 + lambda * eye(n)); %stabilization --> this turns a Gauß-Newton algorithm into a LM-algorithm
        tmp = 1 - x./f_a;
        beta = - transpose(transpose(tmp(:)) * J);
        delta = alpha \ beta;

        %positivity constraints for z-position, signal and background (these must be
        %contained in parameter vector at positions 3, 4 and 5!)
        tmp = a + delta;
        if tmp(3) < 0
            delta(3) = 0;
        end
        if tmp(4) < 0 
            delta(4) = 0; 
        end
        if tmp(5) < 0 
            delta(5) = 0; 
        end

        f_a_new = f(a + delta); 
           
        iter = iter + 1;

        if norm(f_a_new - x) <= norm(f_a - x) 
            lambda = lambda/factor;
            a = a + delta; %parameter update
        else
            lambda = lambda*factor;
        end 

       
        if (norm(f_a_new - f_a) < tol) || (iter >= maxiter)
            break
        end

    end
    %disp(iter)
end