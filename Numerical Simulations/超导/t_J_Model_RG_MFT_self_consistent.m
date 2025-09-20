%[text] [t-J模型的重整化平均场自洽方程的求解（附Julia代码） - 知乎](https://zhuanlan.zhihu.com/p/10019378792)
%[text] 

function [deltas, Delta0s, DeltaSC] = main()
    % 主函数
    N0 = 50;    
    deltas = zeros(N0, 1);
    Delta0s = zeros(N0, 1);
    DeltaSC = zeros(N0, 1);
    
    mu_min = 0.01;
    mu_max = 1.5;
    
    for i = 1:N0
        mu = -(i*(mu_max - mu_min)/N0 + mu_min);
        [delta, Delta0] = RMFT(mu);
        deltas(i) = delta;
        Delta0s(i) = Delta0;
        DeltaSC(i) = g_t(delta)*Delta0;
    end
    
    % 绘图
    figure;
    plot(deltas, [Delta0s, DeltaSC], 'LineWidth', 2);
    xlabel('hole concentration δ');
    ylabel('paring potential Δ');
    legend('Δ0', 'ΔSC');
    grid on;
    set(gcf, 'Position', [100, 100, 300, 400]);
end

function [delta, Delta0] = RMFT(mu)
    % RMFT函数
    t = 5;
    J = 1;
    
    L = 84;
    K = 2*pi/L;
    
    % 初始猜测
    delta = rand();
    Delta0 = rand();
    chi0 = rand();
    
    for step = 1:10000
        gt = g_t(delta);
        gs = g_s(delta);
        
        Delta0_new = 0;
        chi0_new = 0;
        delta_new = 0;
        
        for kx = -pi:K:pi
            for ky = -pi:K:pi
                deltak = delta_k(kx, ky, gs, J, Delta0);
                zetak = zeta_k(kx, ky, gt, gs, t, J, chi0, mu);
                Ek = E_k(zetak, deltak);
                
                Delta0_new = Delta0_new + 1/L^2 * (cos(kx)-cos(ky)) * deltak/(2*Ek);
                chi0_new = chi0_new - 1/L^2 * (cos(kx)+cos(ky)) * zetak/(2*Ek);
                delta_new = delta_new + 1/L^2 * zetak/Ek;
            end
        end
        
        % 自洽精度为delta_new - delta < 10^-12
        if abs(delta_new - delta) < 1e-12
            break;
        end
        
        Delta0 = Delta0_new; 
        chi0 = chi0_new;
        delta = delta_new; 
    end
    
    fprintf('%f\n', delta);
    fprintf('%f\n', Delta0);
end

function res = g_s(delta)
    res = 4/(1+delta)^2;
end

function res = g_t(delta)
    res = (2*delta)/(1+delta);
end

function res = delta_k(kx, ky, gs, J, delta0)
    res = 3*gs*J*delta0/4 * (cos(kx)-cos(ky));
end

function res = zeta_k(kx, ky, gt, gs, t, J, chi0, mu)
    res = -(2*gt*t + 3/4*gs*J*chi0) * (cos(kx)+cos(ky)) - mu;
end

function res = E_k(zeta_k, delta_k)
    res = sqrt(zeta_k^2 + delta_k^2);
end

% 运行主函数并计时
tic;
[deltas, Delta0s, DeltaSC] = main();
toc;


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
