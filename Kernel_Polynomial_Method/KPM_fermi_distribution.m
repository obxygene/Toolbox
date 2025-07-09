function occupation = KPM_fermi_distribution(energy, mu, temperature)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
    if temperature < 0
        error("temperature must be non-negative");
    elseif temperature == 0
        occupation = double(energy < mu);
    else
        occupation = 1./(1 + exp((energy - mu)/temperature)); 
    end
end