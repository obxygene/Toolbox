function Transmission = GreenFunction_Transmission(Gamma_Left, Gamma_Right, G_Centre_Ret)
%
% 利用Fisher-Lee关系通过格林函数计算透射率
% T = Tr[Γ_L * G_C^ret * Γ_R * G_C^adv]
% G_C^ret = G_C^adv'
%   由于输入的格林函数为复数，在使用该函数时一般需要在外面取实部
    Transmission = trace(Gamma_Left * G_Centre_Ret * Gamma_Right * G_Centre_Ret');
end

