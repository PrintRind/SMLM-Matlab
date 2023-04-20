function [a, iter, LLR] = LM_poisson(f, Jac, x, v0, tol, maxiter, lambda, factor)
    %see supplementary materials of 
    %[1] Laurence, Ted A., and Brett A. Chromy. 
    % "Efficient maximum likelihood estimator fitting of histograms." 
    % Nature methods 7.5 (2010): 338-339.
    %f...model function of the molecule, takes the parameters to estimate
    %as inputs
    %x...measured image of the molecule (2D array)

    a = v0'; %initial parameters; IMPORTANT: a must be a column vector!
    n = length(a);
    iter = 0;

    df = numel(x) - numel(a); %number of degrees of freedeom for the Chi² distribution
    %df = 17*17 - 5; %interesting!

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
        alpha = transpose(J) * (J .* repmat(tmp,1,5)); %matches Eq. 9 (left) in suppl. of Ref.[1]
        alpha = alpha .* (1 + lambda * eye(n)); %stabilization --> this turns a Gauß-Newton algorithm into a LM-algorithm
        tmp = 1 - x./f_a;
        beta = - transpose(transpose(tmp(:)) * J); %matches Eq. 9 (right) in [1]
        delta = alpha \ beta; %solving Eq. 10 in suppl. of Ref [1]

        %positivity constraints for z-position, signal and background (these must be
        %contained in parameter vector at positions 3, 4 and 5!)
        tmp = a + delta;
        %if tmp(3) < 0
        %    delta(3) = 0; %positivity constraint for z-position
        %end
        if tmp(4) < 0 
            delta(4) = 0; %positivity constraint for signal
        end
        if tmp(5) < 0 
            delta(5) = 0; %positivity constraint for background
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

    %calculate goodness of fit (log-likelihood ratio) for a later filtering step
    LLR = 2*sum(sum((f_a_new - x .* log(f_a_new)) - (x - x .* log(x))));
    %LLR(isnan(LLR)) = inf; 
    %P = round(chi2cdf(real(LLR),df,'upper')*100); %fraction of molecules (in percent) that will be rejected althogh they fulfill the null-hypothesis 

    %     figure(5);
%     subplot(1,2,1); 
%     imagesc(x);
%     subplot(1,2,2);
%     imagesc(f_a_new);
%     pause;
   

    if iter == maxiter
        disp("max. iterations reached! Fit probably erroneous.")
    end
end