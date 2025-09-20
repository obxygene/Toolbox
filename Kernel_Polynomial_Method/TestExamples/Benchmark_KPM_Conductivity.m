clear
col_tot = 200;
row_tot = 200;
t = 1;
t_x = t;
t_y = t;
epsilon = 0;
phi = 0:0.01:0.2;
Temperature = 0;
    
N_randvec = 1;
N_moments = 500;
N_points = 2*N_moments;

Conductivity_xx = zeros(length(phi), 1);
Conductivity_xy = zeros(length(phi), 1);

% Conductivity_xx = zeros(length(phi), N_points);
% Conductivity_xy = zeros(length(phi), N_points);
    tic;
for jj = 1:length(phi)
    H_sparse = get_Nearest_Square_Hopping_Sparse_V2(row_tot, col_tot, t_x, t_y,[],phi(jj));
    % H_sparse = Hamiltonian_2D_KPM(col_tot, row_tot, epsilon, phi, t);
    x_operator = kron(spdiags((1:col_tot)', 0, col_tot, col_tot), speye(row_tot));
    y_operator = kron(speye(col_tot), spdiags((1:row_tot)', 0, row_tot, row_tot));
    v_operator_x = 1i * sparse(H_sparse * x_operator - x_operator * H_sparse);
    v_operator_y = 1i * sparse(H_sparse * y_operator - y_operator * H_sparse);
    clear x_operator y_operator

%%
    %% KPM计算电导率

    % N_points = 1024;
    [H_tilde, a, b] = KPM_scaleHamiltonian(H_sparse);
    e_scaled = KPM_Chebyshev_abscissas(N_points);
    energies = e_scaled * a + b;
    % mu = linspace(-1,1,101) * a + b;
    mu = -3.7;
%%
    % % 制造[1,0,0...0]向量与kwant对比
    % rand_ket = sparse(zeros(size(H_tilde,1),1));
    % rand_ket(1) = 1;
    % Moments_Mat = KPM_Moments_Correlator(H_tilde, 1, N_moments, v_operator_x, v_operator_x, 'gpu', rand_ket);
    
    Moments_Mat = KPM_Moments_Correlator(H_tilde, N_randvec, N_moments, v_operator_x, v_operator_x, 'gpu');    
    Moments_Kernel = KPM_Kernel_Correction(Moments_Mat, N_moments);   
    integral_factor = zeros(1, length(energies));
    integral = zeros(1, length(mu));
    
    for ii = 1:length(e_scaled)
        Gamma = KPM_Correlator_basis(e_scaled(ii), N_moments);
        integral_factor(ii) = tensorprod(Moments_Kernel, Gamma.','all');
    end
    
    for ii = 1:length(mu)
        fermi_distribution = KPM_fermi_distribution(energies, mu(ii), Temperature);
        integrand = (fermi_distribution ./ ((1 - e_scaled.^2).^2)) .* integral_factor;
        integral(ii) = trapz(e_scaled, integrand);
    end
    % Original coefficient
    % Conductivity_xx = 2 * integral / a^2 / pi^2 ;
    % kwant coefficient
    temp = (8 / a^2) * integral;
    area = 1;
    Conductivity_xx(jj,:) = temp / area;
    clear Moments_Kernel Moments_Mat
%%
    % % 制造[1,0,0...0]向量与kwant对比
    % rand_ket = sparse(zeros(size(H_tilde,1),1));
    % rand_ket(1) = 1;
    % Moments_Mat = KPM_Moments_Correlator(H_tilde, 1, N_moments, v_operator_x, v_operator_x, 'gpu', rand_ket);
    
    Moments_Mat = KPM_Moments_Correlator(H_tilde, N_randvec, N_moments, v_operator_x, v_operator_x, 'gpu');    
    Moments_Kernel = KPM_Kernel_Correction(Moments_Mat, N_moments);    
    integral_factor = zeros(1, length(energies));    
    integral = zeros(1, length(mu));    
    
    for ii = 1:length(e_scaled)
        Gamma = KPM_Correlator_basis(e_scaled(ii), N_moments);
        integral_factor(ii) = tensorprod(Moments_Kernel, Gamma.','all');
    end
    
    for ii = 1:length(mu)
        fermi_distribution = KPM_fermi_distribution(energies, mu(ii), Temperature);
        integrand = (fermi_distribution ./ ((1 - e_scaled.^2).^2)) .* integral_factor;
        integral(ii) = trapz(e_scaled, integrand);
    end
    % Conductivity = 4 * integral / a^2 ;
    Conductivity_xy(jj,:) = 2 * 4 / a^2 * integral ;
end
toc;
%%
% plot(mu, real(Conductivity_xx));
% plot(mu, real(Conductivity_xy));
% plot(phi,real(Conductivity_xy))
% plot(mu, real(Conductivity_xx), mu, real(Conductivity_xy));
plot(phi, real(Conductivity_xx), phi, real(Conductivity_xy));
legend({'$\sigma_{xx}$','$\sigma_{xy}$'},'Interpreter','latex')
function H = Hamiltonian_2D_KPM(Length, Width, epsilon, phi, t)
    
    H_00 = sparse(epsilon * eye(Width) - t * diag(ones(1, Width-1),1) - t * diag(ones(1, Width-1),-1));
    H_01 = sparse(-t * diag(exp(1i * 2*pi* (1:Width) * phi)));

    H = kron(speye(Length), H_00) + kron(spdiags(ones(Length, 1),1, Length, Length), H_01') + kron(spdiags(ones(Length, 1),-1, Length, Length), H_01);
    

end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
