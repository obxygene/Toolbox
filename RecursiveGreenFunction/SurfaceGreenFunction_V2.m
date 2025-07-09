
%% 计算表面格林函数- Lopez Sancho -J. Phys. F: Met. Phys. 15 851

%% 参数介绍
% 输入参数
% H_00 主层层内哈密顿量
% V    主层层间跃迁哈密顿量
% omega 频率ω
% Accuracy 以精度为目标
% varagin{1} 无穷小虚部eta取值，若无输入，则默认为10^(-5)

% 输出参数
% G_00 在频率ω下对应于H00和V的表面格林函数
function [G_00, varargout] = SurfaceGreenFunction_V2(H00, V, omega, Accuracy, varargin)
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
    G_00 = complex(zeros(Width, Width));
    Count = 0;
    Err = Inf;
%%  迭代计算自能
% 初始化自能参数
    alpha = complex(V);
    beta = complex(V');
    Sigma = complex(H00);
    Sigma_tilde = complex(H00);   
    % 迭代
    while abs(Err) > Accuracy
        % Disorder case
        % if nargin > 5
        %     if varargin{2} == "onsite_disorder"
        %         g_temp = ((omega+1i*eta)*Iden - diag((rand(1,Width)-1/2)*varargin{3}) -  Sigma_tilde)\Iden;
        %     end
        % else
        Err_Norm = norm(Sigma,'fro');
        g_temp = ((omega+1i*eta)*Iden - Sigma_tilde)\Iden;
        % end
        Sigma = Sigma + alpha * g_temp * beta;
        Sigma_tilde = Sigma_tilde + beta*g_temp*alpha + alpha*g_temp*beta;
        alpha = alpha * g_temp * alpha;
        beta = beta * g_temp * beta;
        Err = Err_Norm - norm(Sigma,'fro');
        Count = Count + 1;
    end
    % 输出表面格林函数
    G_00 = ((omega+1i*eta)*Iden - Sigma)\Iden;    
    varargout{1} = Count;
end









