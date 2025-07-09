function Moments = KPM_Moments_Operator(H_Scaled, N_randvec, N_moments, Operator, varargin)
%% Kernel Polynomial Method - Moment Calculation with operator
% Calculate the Moments for the expansion coefficient of a function
% input
% H_Scaled[Matrix]: Hamiltonian after scaled to [-1,1]
% N_randvec[Integer]: random vector used to calculate stochastic trace
%                   Generally around 20
% N_moments[Integer]: number of total moments needs to calculate
%                   The number of moments are N_moments in total, but the
%                   first moment is 1, so in this program only last
%                   (N_moments-1) moments are calculated
% varargin{1}[String]: if set as "gpu", the program will calculate using
%                   gpu
% 
% This function used GPU accelerating(Seems not so good)
% 

%% GPU accelerating for large matrix (Seems not so good)
    if nargin > 4 
        if size(H_Scaled,1) > 500
            if canUseGPU
                gpuflag = "gpuArray";
                sprintf("The computer has gpu, for matrix larger than 500 the gpu is initiated")
            else 
                gpuflag = "double";
            end
        else
            gpuflag = "double";
        end
    else
        gpuflag = "double";
    end
%% Main part

    Moments = zeros(1, N_moments, gpuflag);
    Moments(1) = 1;
    for jj = 1:N_randvec
        randvec = exp(2*pi*1i * rand(size(H_Scaled,1), 1,gpuflag));
        randvec = randvec / norm(randvec);
        alpha_minus = randvec;
        alpha = H_Scaled * alpha_minus;
        beta = randvec' * Operator;
        %% Normal calculation
        for kk = 2: N_moments
            Moments(kk) = Moments(kk) + real(beta * alpha);
            alpha_plus = 2 * H_Scaled * alpha - alpha_minus;
            alpha_minus = alpha;
            alpha = alpha_plus;
        end
    end

    % Random average
    Moments(2:end) = real(Moments(2:end)) / N_randvec;

    % If calculation is carried on GPU, gather the data onto memory
    if ~isempty(gpuflag)
        Moments = gather(Moments);
    end
end


%% Old version
% function Moments = KPM_Moments_Calculation(H_Scaled, N_randvec, N_moments, varargin)
% %% Kernel Polynomial Method - Moment Calculation
% % Calculate the Moments for the expansion coefficient of a function
% % input
% % H_Scaled[Matrix]: Hamiltonian after scaled to [-1,1]
% % N_randvec[Integer]: random vector used to calculate stochastic trace
% %                   Generally around 20
% % N_moments[Integer]: number of total moments needs to calculate
% %                   The number of moments are N_moments in total, but the
% %                   first moment is 1, so in this program only last
% %                   (N_moments-1) moments are calculated
% % 
% % This function used GPU accelerating
% % 
% 
% % GPU accelerating for large matrix
%     if size(H_Scaled,1) > 500
%         if canUseGPU
%             gpuflag = "gpuArray";
%             sprintf("The computer has gpu, for matrix larger than 500 the gpu is initiated")
%         else 
%             gpuflag = "double";
%         end
%     else
%         gpuflag = "double";
%     end
% 
%     Moments = zeros(1, N_moments, gpuflag);
% 
% 
%     for jj = 1:N_randvec
%     randvec = exp(2*pi*1i * rand(size(H_Scaled,1), 1,gpuflag));
%     randvec = randvec / norm(randvec);
%     alpha_minus = randvec;
%     alpha = H_Scaled * alpha_minus;
%         for kk = 1: (N_moments - 1)
%             % randvec = exp(2*pi*1i * linspace(0,1,size(H_tilde,1)))*exp(1i * jj);
%             Moments(kk) = Moments(kk) + real(randvec' * alpha);
%             alpha_plus = 2 * H_Scaled * alpha - alpha_minus;
%             alpha_minus = alpha;
%             alpha = alpha_plus;
%         end
%     end
%     Moments = Moments / N_randvec;
%     % If calculation is carried on GPU, gather the data onto memory
%     if ~isempty(gpuflag)
%         Moments = gather(Moments);
%     end
% end

