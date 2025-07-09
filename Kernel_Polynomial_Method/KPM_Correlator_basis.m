function Gamma = KPM_Correlator_basis(energy, N_moments)
%% Kernel Polynomial Method - Basis of Kubo-Bastin formula
% 
%  Γ_{mn}(E)
% = (E - i n \sqrt{1 - E^2}) e^{i n \arccos(E)} T_m(E)
% + (E + i m \sqrt{1 - E^2}) e^{-i m \arccos(E)} T_n(E),
% 
% Notice if chebyshev abscissas was chosen, the acos can be omitted.
    n_array = 0: N_moments - 1;

    T_m = cos(n_array * acos(energy));
    % Not conjugate transpose
    Gamma = T_m.' * ((energy - 1i * n_array * sqrt(1 - energy.^2)) .* exp(1i * n_array * acos(energy)));
    % (m, n) are symmetric
    Gamma = Gamma + Gamma';
end