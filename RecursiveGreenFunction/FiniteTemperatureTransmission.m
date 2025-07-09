
%% 计算有限温度效应的透射率，
% 输入参数
% omega [1*N]: 能量，以k_B为单位
% Transmission_0T [1*N] : 对应能量的零温透射率
% mu [1*1]: 化学势
% Tem [1*1]: 当前温度

function Trans_Tem = FiniteTemperatureTransmission(omega, Transmission_0T, mu, Tem)
% Ensure the dimension of transmission and omega are matched
    if sum(size(omega) ~= size(Transmission_0T))
        Transmission_0T = Transmission_0T';
    end
        Trans_Tem = sum(Transmission_0T .*arrayfun(@(E) fermi_derivative(E-mu,Tem),omega) *(omega(2)-omega(1)));

end

% 局部函数，定义费米分布的导数
function f_derivative = fermi_derivative(E,T)

    f_derivative = 1./(4*T).*sech(E./(2*T)).^2;
end
    