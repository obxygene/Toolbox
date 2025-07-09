function omega_dct = KPM_Chebyshev_abscissas(N_points)
%% Kernel Polynomial Method - abscissas
% give the variable used after discrete cosine transformed function
% 给出chebyshev多项式的节点
%   此处显示详细说明
    omega_dct = cos(pi * ((N_points-1:-1:0)+0.5)/N_points);
end