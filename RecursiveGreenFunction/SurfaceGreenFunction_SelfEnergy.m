function Sigma = SurfaceGreenFunction_SelfEnergy(G_00, H01)
%% SURFACEGREENFUNCTION_SELFENERGY 
% 计算经过表面格林函数修正的自能，输入表面格林函数与边界相互作用，输出自能
% H01定义为中心区与导线区的lead连接处的跃迁哈密顿量
% Σ = H01 * G_00 * H01^\dagger
    Sigma = H01 * G_00 * H01';
end

