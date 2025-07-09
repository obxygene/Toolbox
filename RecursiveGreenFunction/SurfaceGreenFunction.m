
%% 计算表面格林函数- Lopez Sancho -J. Phys. F: Met. Phys. 15 851

%% 参数介绍
% 输入参数
% H_00 主层层内哈密顿量
% V    主层层间跃迁哈密顿量
% omega 频率ω
% Iteration 表面格林函数迭代次数
% varagin{1} 无穷小虚部eta取值，若无输入，则默认为10^(-5)

% 输出参数
% G_00 在频率ω下对应于H00和V的表面格林函数
function G_00 = SurfaceGreenFunction(H00, V, omega, Iteration, varargin)
    %% Parameter setup 参数设定
    % 若输入参数含有自设定无穷小虚部，则采用所设eta，否则默认10^(-4)
    if rem(nargin,2) ~= 0
        eta = varargin{1};
    else
        % 1e-4为测试后相对精确值 
        eta = 1e-4;
    end
    % 通过输入哈密顿量的行数判断单层原子有多少个，定义为Width
    Width = size(H00,1);
    Iden = eye(Width);
    G_00 = zeros(Width, Width);
    
%%  迭代计算自能
% 初始化自能参数
    alpha = V;
    beta = V';
    Sigma = H00;
    Sigma_tilde = H00;   
    % 迭代
    for jj = 1:Iteration
        if nargin > 5
            if varargin{2} == "onsite_disorder"
                g_temp = ((omega+1i*eta)*Iden - diag((rand(1,Width)-1/2)*varargin{3}) -  Sigma_tilde)\Iden;
            end
        else
                g_temp = ((omega+1i*eta)*Iden - Sigma_tilde)\Iden;
        end
        Sigma = Sigma + alpha * g_temp * beta;
        Sigma_tilde = Sigma_tilde + beta*g_temp*alpha + alpha*g_temp*beta;
        alpha = alpha * g_temp * alpha;
        beta = beta * g_temp * beta;
    end
    % 输出表面格林函数
    G_00(:,:) = ((omega+1i*eta)*Iden - Sigma)\Iden;
    
    %% 辅助debug部分，若无必要可注释
    % tfin = toc(tstart);
    % figure(2)
    % plot(reshape(-imag(G_00(1,1,:)/pi),[1,N_point]),omega)
    % fprintf('The time used in surface green function for frequency %4.2f is %4.2f second\n',omega, tfin)

end









