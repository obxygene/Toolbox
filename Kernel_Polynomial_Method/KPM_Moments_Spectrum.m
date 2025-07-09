function Moments = KPM_Moments_Spectrum(H_Scaled, N_randvec, N_moments, varargin)
%% Kernel Polynomial Method - Moment Calculation
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
    if nargin > 3 
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

    for jj = 1:N_randvec
        randvec = exp(2*pi*1i * rand(size(H_Scaled,1), 1,gpuflag));
        randvec = randvec/norm(randvec);
        alpha_minus = randvec;
        alpha = H_Scaled * alpha_minus;

        %% Normal calculation
        % Moments(1) = 1;
        % for kk = 2: N_moments
        %     Moments(kk) = Moments(kk) + real(randvec' * alpha);
        %     alpha_plus = 2 * H_Scaled * alpha - alpha_minus;
        %     alpha_minus = alpha;
        %     alpha = alpha_plus;
        % end

        %% Fast Calculation
        Moments_temp = zeros(1, N_moments, gpuflag);
        Moments_temp(1) = randvec' * randvec;
        Moments_temp(2) = randvec' * alpha; 
        for kk = 2: (N_moments/2)
            alpha_plus = 2 * H_Scaled * alpha - alpha_minus;
            Moments_temp(2 * kk - 1) = 2 * (alpha' * alpha) - Moments_temp(1);
            Moments_temp(2 * kk) = 2 * (alpha_plus' * alpha) - Moments_temp(2);
            alpha_minus = alpha;
            alpha = alpha_plus;
        end
        if mod(N_moments, 2)
            Moments_temp(end) = 2 * (alpha_plus' * alpha_plus) - Moments_temp(1);
        end

        % prepare to do the average
        Moments = Moments + Moments_temp;
    end

    % Random average
    Moments = real(Moments)/N_randvec;

    % If calculation is carried on GPU, gather the data onto memory
    if gpuflag == "gpuArray"
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

