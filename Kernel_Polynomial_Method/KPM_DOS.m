function [omega, DOS] = KPM_DOS(H_sparse, varargin)
%% Kernel Polynomial Method - DOS calculation
% 
% input:
% H_sparse[Matrix]: input sparse Hamiltonian
% 
% default epsilon = 0.05 
% default ScaleFactor = 1
    if nargin > 1
        N_randvec = varargin{1};
        if nargin > 2
            N_moments = varargin{2};
            if nargin > 3
                N_points = varargin{3};
            else
                N_points = 2 * N_moments;
            end
        else
            N_moments = 200;
        end
    else
        N_randvec = 20;
        N_moments = 200;
        N_points = 2 * N_moments;
    end
    [H_tilde, a, b] = KPM_scaleHamiltonian(H_sparse);
    Moments = KPM_Moments_Spectrum(H_tilde, N_randvec, N_moments);
    Moments_Kernel = KPM_Kernel_Correction(Moments, N_moments);
    [omega_dct, DOS] = KPM_DCT(Moments_Kernel, N_points);
    DOS = sqrt(length(Moments_Kernel)) * DOS;
    omega = KPM_rescale(omega_dct, a, b);

end










