
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



kx_step = 100;
kx_list = linspace(0, 2*pi/sqrt(3), kx_step);
t_1 = 1;
t_2 = 1/3;
phi = -pi/4;
M = 2/3;
eigen_list = zeros(kx_step,y_num);
for i = 1:kx_step
    temp = Haldane_zigzag_Edge_Hamiltonian(kx_list(i), ny_unit_cell, t_1, t_2, phi, M);
    eigen_list(i,:) = eig(temp);
end
figure(1)
plot(kx_list,eigen_list)




%% Edge state at kx = 0

    EdgeHamiltonian = Haldane_zigzag_Edge_Hamiltonian(2, ny_unit_cell, t_1, t_2, phi, M);
    [Wave, ~] = eig(EdgeHamiltonian);
    rho = abs(Wave(1:4:y_num,:)).^2 + abs(Wave(4:4:y_num,:)).^2 + ...
            abs(Wave(2:4:y_num,:)).^2+abs(Wave(3:4:y_num,:)).^2;
    figure(2)
    plot(rho);




%% function definition
function H = Haldane_zigzag_Edge_Hamiltonian(kx, y_num, t_1, t_2, phi, M)
% Haldane哈密顿量，为得到x方向的边缘态，切开x方向后的哈密顿量
% 参考资料来自于https://phys-drp.mit.edu/sites/default/files/thaodinh.pdf
% kx, ky为布里渊区内的晶格动量,t为跃迁强度
% 在这个程序中，为了得到边缘态，假设在y方向上有穷而x方向上无穷，则
    k = [kx, 0];
    a_1 = [0,1]';
    a_2 = [-sqrt(3)/2, -1/2]';
    a_3 = [sqrt(3)/2 , -1/2]';
    % a = [a_1, a_2, a_3];
    % x_trans = [3,0]';
    b_1 = [sqrt(3), 0]';
    b_2 = [-sqrt(3)/2, -3/2]';
    b_3 = [-sqrt(3)/2, 3/2]';
    % b = [b_1, b_2, b_3];

    H_0 =        [M, t_1, t_2 * exp(1i * phi), 0;
                  t_1,-M, t_1, t_2 * exp(1i * phi);
                  t_2 * exp(-1i * phi), t_1, M, t_1;
                  0, t_2 * exp(-1i * phi), t_1, -M];

    H_12 =       [0, 0, 0, 0;
                  0, 0, 0, 0;
                  t_2 * exp(-1i * phi), 0, 0, 0;
                  t_1, t_2 * exp(-1i * phi), 0, 0];

    H_21 =       H_12';

    T_0 =  zeros(4,4);
    T_0(1,1) = t_2 * exp(-1i * phi);
    
    T_0(2,1) = t_1; 
    T_0(2,2) = t_2 * exp(+1i * phi);
    T_0(2,4) = t_2 * exp(-1i * phi);
    
    T_0(3,1) = t_2 * exp(1i * phi);
    T_0(3,3) = t_2 * exp(-1i * phi);
    T_0(3,4) = t_1;

    T_0(4,4) = t_2 * exp(1i * phi);
    
    TR_12 =      [0, 0, 0, 0;
                 0, 0, 0, 0;
                 t_2 * exp(1i * phi), 0, 0, 0;
                 0, 0, 0, 0];

    TR_21 =      [0, 0, 0, 0;
                 0, 0, 0, t_2*exp(-1i*phi);
                 0, 0, 0, 0;
                 0, 0,0, 0];

    TL_12 =      [0, 0, 0, 0;
                 0, 0, 0, 0;
                 0, 0, 0, 0;
                 0, t_2 * exp(1i * phi), 0, 0];

    TL_21 =      [0, 0, t_2 * exp(-1i*phi), 0;
                 0, 0, 0, 0;
                 0, 0, 0, 0;
                 0, 0,0, 0];

    
    % Onsite part
    H =   kron(diag(ones(y_num,1)), H_0) ...
        + kron(diag(ones(y_num - 1,1),1), H_12) ...
        + kron(diag(ones(y_num - 1,1),-1), H_21);
    
    % Right hopping
    H = H + (kron(diag(ones(y_num,1)), T_0) ...
        + kron(diag(ones(y_num - 1,1),1), TR_12) ...
        + kron(diag(ones(y_num - 1,1),-1), TR_21)) * exp(1i * k * b_1);

    % Left hopping
    H = H + (kron(diag(ones(y_num,1)), T_0') ...
        + kron(diag(ones(y_num - 1,1),1), TL_12) ...
        + kron(diag(ones(y_num - 1,1),-1), TL_21)) * exp(-1i * k * b_1);

    % Boundary hopping 


end













