function [omega_dct, Reconstruct] = KPM_DCT(Moments_Kernel, N_points,varargin)
%% Kernel Polynomial Method - Discrete Cosine Transform
% 利用chebyshev多项式的节点，采用离散余弦变换的方法快速计算重构的函数

    omega_dct = KPM_Chebyshev_abscissas(N_points);
    % matlab的离散余弦变换dct定义中，第一项系数与常用的相差根号2
    Moments_Kernel(1) = Moments_Kernel(1)/sqrt(2);
    % 若对象是矩阵，则作两次dct变换
    if ~isvector(Moments_Kernel)
        Reconstruct = flip(dct(Moments_Kernel,N_points,1,'Type', 3));
        Reconstruct = flip(dct(Reconstruct,N_points,2,'Type', 3));
    else
        % 若对象是向量，只作一次dct变换
        Reconstruct = flip(dct(Moments_Kernel,N_points,'Type', 3))./(pi * sqrt(1 - omega_dct.^2));
    end
end