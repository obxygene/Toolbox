function KernelCoeff = KPM_Kernel(N_moments, Kernel_Type, varargin)
%% Kernel Polynomial Method - Kernel coefficient
% generate the coefficient that can suppress the gibbs phenomenon


if Kernel_Type == "Jackson"
    KernelCoeff = ((N_moments+1-(0:N_moments-1)).*cos(pi * (0:N_moments-1)/(N_moments+1)) + ...
                   sin(pi * (0:N_moments-1)/(N_moments+1))*cot(pi/(N_moments+1)) )/(N_moments + 1);
elseif Kernel_Type == "Lorentz"
    if nargin < 3
        lambda = 3;
    else
        lambda = varargin{1};
    end
    KernelCoeff = sinh(lambda*(1-(0:N_moments-1)/N_moments))/sinh(lambda);
else
    error("Kernel_Type should be either 'Jackson' or 'Lorentz'");
end