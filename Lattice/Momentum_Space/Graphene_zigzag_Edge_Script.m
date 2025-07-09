
clear;
%% Parameter Setup
% y_num is the total number of atoms in y direction
% N_unitcell is the number of atoms in one unit cell, zig zag direction
%   have 4 atoms in one unit cell
% unitcell_atom_y is the number of y directional atoms in one unit cell.
%   which is different from the N_unitcell because in armchair direction
%   this quantity will be 2

% 注意到仅在x方向上应用Bloch定理时，每个单胞之间的间隔是sqrt(3)，
% 故布里渊区应当是2pi/sqrt(3)

y_num = 100;
N_unitcell = 4;
unitcell_atom_y = 4;
ny_unit_cell = y_num / unitcell_atom_y;


kx_step = 1000;
kx_list = linspace(-2*pi/sqrt(3), 2*pi/sqrt(3), kx_step);
t = 1;
M = 0 * t;
eigen_list = zeros(kx_step,y_num);
for i = 1:kx_step
    temp = Graphene_zigzag_Edge_Hamiltonian(kx_list(i), ny_unit_cell, t, M);
    eigen_list(i,:) = eig(temp);
end
figure
plot(kx_list,eigen_list)




%% Edge state at kx = 0
figure
    EdgeHamiltonian = Graphene_zigzag_Edge_Hamiltonian(3.5, ny_unit_cell, t, M);
    [Wave, ~] = eig(EdgeHamiltonian);
    rho = abs(Wave(1:4:y_num,:)).^2 + abs(Wave(4:4:y_num,:)).^2 + ...
            abs(Wave(2:4:y_num,:)).^2+abs(Wave(3:4:y_num,:)).^2;
    figure(2)
    plot(rho);


%% function definition
function H = Graphene_zigzag_Edge_Hamiltonian(kx, y_num, t, M)
% 石墨烯哈密顿量，为得到x方向的边缘态，切开x方向后的哈密顿量
% 参考资料来自于https://phys-drp.mit.edu/sites/default/files/thaodinh.pdf
% kx, ky为布里渊区内的晶格动量,t为跃迁强度
% 在这个程序中，为了得到边缘态，假设在y方向上有穷而x方向上无穷，则
    k = [kx, 0];
    a_1 = [0,1]';
    a_2 = [-sqrt(3)/2, -1/2]';
    a_3 = [sqrt(3)/2 , -1/2]';
    % a = [a_1, a_2, a_3];
    % x_trans = [3,0]';
    % b_1 = [0, sqrt(3)]';
    % b_2 = [-3/2, -sqrt(3)/2]';
    % b_3 = [+3/2, -sqrt(3)/2]';
    % b = [b_1, b_2, b_3];

    H_0 =        [M, t, 0, 0;
                  t,-M, t, 0;
                  0, t, M, t;
                  0, 0, t, -M];

    H_12 =       [ 0, 0, 0, 0;
                   0, 0, 0, 0;
                   0, 0, 0, 0;
                   t, 0, 0, 0];

    H_21 =       [ 0, 0, 0, t;
                   0, 0, 0, 0;
                   0, 0, 0, 0;
                   0, 0, 0, 0];

    T_0 =        [0, t * exp(1i * k * a_2), 0, 0;
                  t * exp(-1i * k * a_2), 0, 0, 0;
                  0, 0, 0, t * exp(1i * k * a_3);
                  0, 0, t * exp(-1i * k * a_3), 0];

    % boundary_block =zeros(N_unitcell, N_unitcell);
    % H_unit = [H_0 + T_0, H_12;
    %           H_21, H_0 + T_0];

    H =   kron(diag(ones(y_num,1)), ((H_0+T_0))) ...
        + kron(diag(ones(y_num - 1,1),1), H_12) ...
        + kron(diag(ones(y_num - 1,1),-1), H_21);

    % Boundary hopping 


end















