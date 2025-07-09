function Moments_Kernel = KPM_Kernel_Correction(moments, N_moments, Kernel_type, varargin)
%% Kernel Polynomial Method - Kernel coefficient suppress Gibbs phenomenon
% input
% moments[vector/matrix]: unmodified chebyshev coefficients
% 
% N_moments[Integer]: number of total moments needs to calculate, including
%                     0th term
% Kernel_type[str]: choose the kernel type e.g. Jackson/Lorentz
%                   default as Jackson
% varargin{1}[scalar]: the lorentz kernel has default kernel coefficient
    if nargin > 2
        if nargin > 3
            Kernel_Coeff = KPM_Kernel(N_moments, Kernel_type, varargin{1});
        else
            Kernel_Coeff = KPM_Kernel(N_moments, Kernel_type);
        end
    else
        Kernel_Coeff = KPM_Kernel(N_moments, 'Jackson');
    end
    if isvector(moments)
        Moments_Kernel = Kernel_Coeff .* moments;
    elseif ismatrix(moments)
        % Without discrete fourier transform, the coefficient must be
        % divided by 2
        Kernel_Coeff(1) = Kernel_Coeff(1) / 2;
        Moments_Kernel = (Kernel_Coeff' .* Kernel_Coeff) .* moments;
    end
end