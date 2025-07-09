function omega = KPM_rescale(omega_dct, a, b)
%% Kernel Polynomial Method - Rescale the spectrum
% Due to the scaling of Hamiltonian, the spectrum is suppressed onto [-1,1]
% To recover the original spectrum, the variable need to be rescaled
% 此处显示详细说明
    omega = omega_dct*a+b;
end