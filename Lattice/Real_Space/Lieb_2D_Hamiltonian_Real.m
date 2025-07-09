
function H = Lieb_2D_Hamiltonian_Real(epsilon,J_plus,J_minus_x,J_minus_y,N_row,N_col, varargin)
%LIEB_HAMILTONIAN_REAL 给出Lieb晶格的实空间紧束缚哈密顿量
%   参数定义 Parameter definition
% epsilon: 1*3 vector, onsite energy, including chemical potential
% n: label of unit cell
% t_minus: A_n->B_n
% t_y: A_n->C_n
% t_plus: B_n->A_{n+1}

%% Lieb Geometry
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...
% .
% .
% .
%% H_intra 原胞内的跃迁
% J_plus: A_n->B_n and A_n->C_n
% J_minus_x: B_xn->A_xn+1
% J_minus_y: C_yn->A_yn+1
% n: label of unit cell
    H_intra = [epsilon(1), J_plus, J_plus;
               J_plus, epsilon(2),0;
               J_plus, 0, epsilon(3)];
    N_Lieb_site = 3;
%% 原胞间跃迁

    % 构造纵向跃迁元 Transition between row
    temp = get_Nearest_Square_Hopping_Sparse(N_col, N_row, J_minus_y, 0);
    temp_vect = [0,0,1;0,0,0;0,0,0];
    temp = tril(kron(full(temp),temp_vect));
    % 增加纵向的周期性边界条件
    if nargin > 6
        temp(sub2ind(size(temp), N_Lieb_site*N_row:N_Lieb_site*N_row:N_Lieb_site*N_row*N_col,1:N_Lieb_site*N_row:(N_Lieb_site*N_row*(N_col-1)+1))) = J_minus_y;
    end
    H_inter_y = temp + temp';

    % 构造横向跃迁元 Transition between column

    temp = get_Nearest_Square_Hopping_Sparse(N_col, N_row, 0, J_minus_x);
    temp_hori = [0,1,0;0,0,0;0,0,0];
    temp = tril(kron(full(temp),temp_hori));
    H_inter_x = temp + temp';
    


    H = kron(eye(N_row*N_col), H_intra) + H_inter_x + H_inter_y;
end
