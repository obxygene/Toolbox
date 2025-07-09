function Gamma = SurfaceGreenFunction_Broadening(Self_energy_matrix)

% 表面格林函数，输入自能矩阵，返回展宽矩阵
%   此处显示详细说明
    Gamma = 1i * (Self_energy_matrix - Self_energy_matrix');
end

