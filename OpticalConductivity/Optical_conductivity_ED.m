function sigma = Optical_conductivity_ED(H, x, omega, mu, eta, Tem)
arguments (Input)
    H (:,:) {mustBeMatrix}
    x {mustBeVector}
    omega double {mustBeVector} = linspace(-2,2,101)
    mu double {mustBeFloat} = 0
    eta double {mustBeFloat} = 1e-3
    Tem double {mustBeFloat} = 0.1
end

arguments (Output)
    sigma {mustBeVector};
end

sigma = zeros(1, length(omega));
Lorentz_delta = @(epsilon) eta.^2 ./(epsilon.^2 + eta^2)/pi;
x = diag(x);
[eigvecs, eigvals] = eig(H);
eigvals = diag(eigvals);
fermi_dist_vals = 1 ./ (exp((eigvals - mu) / Tem) + 1);
matrix_elements = abs(eigvecs' * x * eigvecs).^2;
for kk = 1:length(omega)
    for ii = 1:length(eigvals)
        E_a = eigvals(ii);
        f_a = fermi_dist_vals(ii);
        omega_ba = eigvals - E_a;
        f_b = fermi_dist_vals;
        delta_vals = Lorentz_delta(omega(kk) - omega_ba);
        sigma(kk) = sigma(kk) + sum(matrix_elements(ii, :)' .* omega_ba .* (f_a - f_b) .* delta_vals);
    end
end    
sigma = sigma / length(x);
end
