function H = HaldaneHamiltonian(kx, ky, t_1, t_2, phi)
% Haldane哈密顿量，通过考虑石墨烯带有相位的次近邻相互作用，从而得到gap
% 参考资料来自于https://phys-drp.mit.edu/sites/default/files/thaodinh.pdf
% kx, ky为布里渊区内的晶格动量，t_1和t_2为跃迁强度，phi为相位
    k = [kx, ky];
    a_1 = [1,0]';
    a_2 = [-1/2, +sqrt(3)/2]';
    a_3 = [-1/2, -sqrt(3)/2]';
    a = [a_1, a_2, a_3];

    b_1 = [0, sqrt(3)]';
    b_2 = [-3/2, -sqrt(3)/2]';
    b_3 = [+3/2, -sqrt(3)/2]';
    b = [b_1, b_2, b_3];
    H11_phase = 0;
    H12_phase = 0;

    H22_phase = 0;
%% Hamiltonian in reference https://phys-drp.mit.edu/sites/default/files/thaodinh.pdf
    % for ii = 1:length(a)
    %     H11_phase = H11_phase + exp(1i * k * b(:,ii));
    %     H12_phase = H12_phase + exp(1i * k * a(:,ii));
    % end

    % H = [t_2 * (exp(1i * phi) * conj(H11_phase) + exp(-1i * phi) * H11_phase), t_1 * conj(H12_phase); ...
    %      t_1 * H12_phase, t_2 * (exp(-1i * phi) * conj(H11_phase) + exp(1i * phi) * H11_phase)];

%% Hamiltonian in 10.1103/PhysRevB.104.045103
    for ii = 1:length(a)
        H11_phase = H11_phase + cos(phi - k * b(:,ii));
        H12_phase = H12_phase + exp(1i * k * a(:,ii));
        H22_phase = H22_phase + cos(phi + k * b(:,ii));
    end



    H = [2*t_2 * H11_phase, t_1 * conj(H12_phase);
         t_1 * H12_phase, 2*t_2 * H22_phase];


end

