function Moments_Mat = KPM_Moments_Correlator(H_Scaled, N_randvec, N_moments, Operator_a, Operator_b, varargin)
%% Kernel Polynomial Method - Correlator calculation using Kubo-Bastin
% Calculate the Correlator

%  χ_{αβ}(µ, T) =
%       \int_{-\infty}^{\infty}dE f(µ-E, T)
%       \left(O_α ρ(E) O_β \frac{dG^{+}}{dE} - O_α frac{dG^{-}}{dE} O_β ρ(E)\right),

% The correlator can be calculated using 2D KPM

% input
% [1] H_Scaled[Matrix]: Hamiltonian after scaled to [-1,1]
% [2] N_randvec[Integer]: random vector used to calculate stochastic trace
%                   Generally around 20
% [3] N_moments[Integer]: number of total moments needs to calculate
% [4] [5]Operator_a/b[Martix]: different operators
%                   if velocity operator, label for different direction a and b
% [6] varargin{1}[String]: if set as "gpu", the program will calculate using
%                   gpu
% [7] varargin{2}[Vector]: if nargin>6, this is the given random vector
% This function used GPU accelerating(Seems not so good)
% 

%% GPU accelerating for large matrix (Seems not so good)
    if nargin > 5
        if size(H_Scaled,1) > 500
            if canUseGPU
                gpuflag = "gpuArray";
                % disp("This computer has gpu, for matrix larger than 500 the gpu is initiated");
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

    Moments_Mat = zeros(N_moments, N_moments, gpuflag);

    for jj = 1:N_randvec
        % 
        if nargin > 6
            rand_ket = varargin{2};
        else
            % Construct Random vectors
            rand_ket = exp(2*pi*1i * randn(size(H_Scaled,1), 1,gpuflag));
            rand_ket = rand_ket / norm(rand_ket);
        end
        rand_ket = rand_ket / norm(rand_ket);
        Kets = zeros(size(H_Scaled,1), N_moments);
        Bras = zeros(size(H_Scaled,1), N_moments);
        %% Calculate T_m(H)\vket{r}
        alpha_minus = rand_ket;
        alpha = H_Scaled * alpha_minus;
        Kets(:, 1) = alpha_minus;
        Kets(:, 2) = alpha;
        for mm = 3: N_moments
            alpha_plus = 2 * H_Scaled * alpha - alpha_minus;
            alpha_minus = alpha;
            alpha = alpha_plus;
            Kets(:, mm) = alpha;
        end
        % Calculate O_b T_n(H)\vket{r}
        Kets = Operator_b * Kets;
        %% Calculate \vbra{r} O_b T_n(H)
        beta_minus = Operator_a * rand_ket;
        beta = H_Scaled * beta_minus;
        Bras(:, 1) = beta_minus;
        Bras(:, 2) = beta;
        for nn = 3: N_moments
            beta_plus = 2 * H_Scaled * beta - beta_minus;
            beta_minus = beta;
            beta = beta_plus;
            Bras(:, nn) = beta;
        end
        Moments_Mat = Moments_Mat + real(Bras' * Kets);
    end

    % Random average
    Moments_Mat = Moments_Mat / N_randvec;
    % If calculation is carried on GPU, gather the data onto memory
    if gpuflag == "gpuArray"
        Moments_Mat = gather(Moments_Mat);
    end
end
