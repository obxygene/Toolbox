

kx_step = 100;
ky_step = 100;
kx_list = linspace(-pi, pi, kx_step);
ky_list = linspace(-pi, pi, ky_step);
t_1 = 1;




energy_band = zeros(kx_step, ky_step, 2);

for i = 1:length(kx_list)
    for j = 1:length(ky_list)
        energy_band(i, j, :) = eig(GrapheneHamiltonian(kx_list(i), ...
            ky_list(j),t_1));
    end
end


surf(kx_list, ky_list, energy_band(:,:,1))
hold on;
surf(kx_list, ky_list, energy_band(:,:,2))


function H = GrapheneHamiltonian(kx, ky, t_1)
%HALDANEHAMILTONIAN 此处显示有关此函数的摘要
%   此处显示详细说明
    k = [kx, ky];
    a_1 = [1,0]';
    a_2 = [-1/2, +sqrt(3)/2]';
    a_3 = [-1/2, -sqrt(3)/2]';
    a = [a_1, a_2, a_3];

    % b_1 = [0, sqrt(3)]';
    % b_2 = [-3/2, -sqrt(3)/2]';
    % b_3 = [+3/2, -sqrt(3)/2]';
    % b = [b_1, b_2, b_3];
    % H11_phase = 0;
    H12_phase = 0;

    for ii = 1:length(a)
        % H11_phase = H11_phase + exp(1i * k * b(:,ii));
        H12_phase = H12_phase + exp(1i * k * a(:,ii));
    end

    H = [0, t_1 * conj(H12_phase); ...
         t_1 * H12_phase, 0];


end