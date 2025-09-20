%[text] [用Matlab求解平均场方程:从自由费米子到BCS平均场理论 - 知乎](https://zhuanlan.zhihu.com/p/1907686051637626370)
%[text] ### 一维费米子
%[text] 
%[text] 考虑一维自由费米子体系
%[text] $\\hat{H}=\\sum\_k \\epsilon\_k \\hat{c}\_k^\\dagger\\hat{c}\_k,\\quad\\epsilon\_k=-2t\\cos k\_x-\\mu$
%[text] 对应的自由能为
%[text] $F=-T\\sum\_k\\ln \\left(1+e^{-\\epsilon\_k/T}\\right)$
%[text] 从而电子数密度由自由能密度对化学势的导数所决定
%[text] $n\_c=-\\frac{1}{N\_s}\\frac{\\partial F}{\\partial \\mu}=\\frac{1}{N\_s}\\sum\_k f\_F(\\epsilon\_k)=\\int^{\\pi}\_{-\\pi}\\frac{dk}{2\\pi}f\_F(\\epsilon\_k),\\quad f\_F(x)=\\frac{1}{e^{x/T}+1}$
%[text] 在绝大多数的情况下，我们更关心固定粒子数密度$n\_c${"editStyle":"visual"}时候的性质，因此我们就需要对于给定的粒子数密度求出对应的化学势，这要求求解上面这个方程。
%[text] 在Matlab中求解一般的方程需要用到fsolve这个命令，其调用格式一般是
%[text] ```matlabCodeExample
%[text] [x, fval] = fsolve('fun', x0, option)
%[text] ```
%[text] 
%[text] x是方程的解，fun是所求方程的名字，X0是方程的初值，option是一些选项，多数时候可以不写。对于我们这里关于nc的方程，设其方程名为F1，那么求解代码就是
%[text] ```matlabCodeExample
%[text] x0=[0.1];
%[text] [x,fval]=fsolve(@F1,x0);
%[text] ```
%[text] 注意方程中还有关于$k\_x${"editStyle":"visual"}的积分，因此还需要用到积分命令quad,从而方程应当改写为
%[text] $\\int\_{-\\pi}^{\\pi}\\frac{dk}{2\\pi}f\_F(\\epsilon\_k)-n\_c=0$
%[text] 从而上述方程对应的函数代码应当为
function [y] = F1(mu, param)
 S = 1 / (2 * pi);
    function [z1] = nf1(kx)
        ek = -2*param.t*cos(kx) - mu;
        z1 = 1 ./ (exp(ek / param.T) + 1);
    end
    Sp1 = integral(@nf1, -pi, pi);
    y = S * Sp1 - param.nc;
end
%[text] 进而可以通过写一个主程序来调用F1，计算给定化学势
param.t = 1;
param.nc = 0.5;
param.T = 0.01 * param.t;
x0 = [0.1];
Fun = @(x)F1(x, param);
[x, fval] = fsolve(Fun, x0) %[output:3b6e727c] %[output:04b8accc] %[output:471c8b91]
%[text] 这其实对应于半满填充的自由费米子模型其化学势$\\mu =0${"editStyle":"visual"}
%[text] ## 二维费米子
%[text] 现在，我们可以考虑二维费米子的情况，为简单期间考虑的是正方晶格体系，这样
%[text] $\\epsilon\_k = -2t(\\cos k\_x+\\cos k\_y)-\\mu$
%[text] 此时平均场方程变为
%[text] $n\_c=-\\frac{1}{N\_s}\\frac{\\partial F}{\\partial \\mu}=\\frac{1}{N\_s}\\sum\_k f\_F(\\epsilon\_k)=\\int^{\\pi}\_{-\\pi}\\frac{dk\_x}{2\\pi} \\int^{\\pi}\_{-\\pi}\\frac{dk\_y}{2\\pi}\\,f\_F(\\epsilon\_k)$
%[text] 给出，现在积分方程需要变为2D情形
function [y] = F2(mu, param)
 S = 1 / (2 * pi)^2;
    function [z1] = nf2(kx, ky)
        ek = -2*param.t*(cos(kx) + cos(ky)) - mu;
        z1 = 1 ./ (exp(ek / param.T) + 1);
    end
    Sp1 = integral2(@nf2, -pi, pi, -pi, pi);
    y = S * Sp1 - param.nc;
end
%[text] 沿用之前的参数，可以求解得到
x0 = [0.1];
Fun = @(x)F2(x, param);
[x, fval] = fsolve(Fun, x0) %[output:5fe99f88] %[output:1a5a592c] %[output:20b16d5a]
%%
%[text] ## 超导BCS平均场方程
%[text] 
%[text] 超导BCS理论在Bogoliubov变换后哈密顿量为
%[text] $\\hat{H}=\\sum\_k E\_k(\\hat{\\alpha}\_k^\\dagger\\hat{\\alpha}\_k+\\hat{\\beta}\_k^\\dagger\\hat{\\beta}\_k)+\\frac{V}{g}|{\\Delta}|^2+\\sum\_k(\\xi\_k-E\_k)$
%[text] 其中准粒子能量为
%[text] $E\_k=\\sqrt{\\xi\_k^2+|\\Delta\_k|^2}$
%[text] 对应的自由能为
%[text] $F=\\frac{V}{g}\\Delta^2+\\sum\_k(\\xi\_k-E\_k)-2T\\sum\_k\\ln\\left(1+e^{-\\beta E\_k}\\right)$
%[text] 带入$\\frac{\\partial E\_k}{\\Delta^2}=\\frac{\\gamma\_k^2}{2E\_k}$
%[text] 可以求得超导能隙方程
%[text] $1=\\frac{g}{V}\\sum\_k\\frac{\\gamma\_k^2}{2E\_k}\\tanh\\frac{E\_k}{2T}$
%[text] 而总电子数密度由如下方程给出
%[text] $n\_c=-\\frac{1}{V}\\frac{\\partial F}{\\partial \\mu}=-\\frac{1}{V}\\left\[\\frac{\\partial \\xi\_k}{\\partial{\\mu}}-\\left(1-2f\_F(E\_k)\\right)\\frac{\\partial E\_k}{\\partial \\mu}\\right\]=\\frac{1}{V}\\sum\_k\\left(1-\\frac{\\xi\_k}{E\_k}\\tanh\\frac{E\_k}{2T}\\right)$
%[text] 所以现在需要求解两个平均场方程，对应的初值也有两个
param.t = 1;
param.nc = 0.5;
param.T = 0.01 * param.t;
param.t1 = 0;
param.t2 = 0;
param.g = 2.0 * t;

% 第一个初值为能隙，第二个初值为化学势
x0 = [0.01, 0];
Fun = @(x)SC_F1(x, param);
[x, ~] = fsolve(Fun, x0) %[output:5b7d4a6e] %[output:1a112624]
%[text] 对应的函数为
function [y] = SC_F1(A, param)
S = 1/(2*pi)^2;
dd = 1e-6;

    function [z1] = Gap_eq(kx, ky)
        xi_k = -2*param.t*(cos(kx)+cos(ky)) - 4*param.t1*cos(kx).*cos(ky)-...
            4*param.t2*(cos(kx).^2+cos(ky).^2-1) - A(2);
        % xi_k = epsilon_k - \mu
        E_k = sqrt(xi_k.^2 + A(1)^2 * gamma(kx,ky).^2+dd);
        % E_k = \sqrt{xi_k^2 + \Delta^2}
        tp1 = (gamma(kx,ky).^2).*tanh(E_k/(2*param.T))./(2*E_k);
        
        z1 = S*tp1;
    end

Sp1 = integral2(@Gap_eq,-pi,pi,-pi,pi);

y1 = param.g * Sp1-1;

    function [z2] = number_density_eq(kx, ky)
        xi_k = -2 * param.t * (cos(kx) + cos(ky)) - 4*param.t1*cos(kx).*cos(ky)-...
            4*param.t2*(cos(kx).^2+cos(ky).^2-1)-A(2);
        E_k = sqrt(xi_k.^2 + A(1)^2 * gamma(kx,ky).^2 + dd);
        
        tp1 = (1- (xi_k./E_k) .*tanh(E_k/(2*param.T)));
        
        z2 = S*tp1;
    end

Sp2 = integral2(@number_density_eq,-pi,pi, -pi,pi);
y2 = Sp2 - param.nc;

y=[y1;y2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = gamma(kx,ky)
%y=cos(kx)-cos(ky); % dx^2-y^2-wave:cos(kx)-cos(ky).   
y=1;               % Uniform s-wave:1.  
%y=cos(kx)+cos(ky); % Extended s-wave:cos(kx)+cos(ky).
%y=sin(kx).*sin(ky);% dxy-wave:sin(kx)*sin(ky).
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:3b6e727c]
%   data: {"dataType":"text","outputData":{"text":"\n<a href = \"matlab: helpview('optim','eqn_solved','CSHelpWindow');\">方程已解<\/a>。\n\nfsolve 已完成，因为按照<a href = \"matlab: helpview('optim','fcn_tolerance_fsolve','CSHelpWindow');\">函数容差<\/a>的值衡量，\n函数值向量接近于零，并且按照梯度的值衡量，\n<a href = \"matlab: helpview('optim','appears_regular','CSHelpWindow');\">问题似乎为正则问题<\/a>。\n\n<<a href = \"matlab: createExitMsg({'optim:fsolve:Exit1basic','fsolve'},{'optim:fsolve:Exit1detailed','5.428370e-14','1.000000e-06','1.163223e-25','1.000000e-03'},true,true);;\">停止条件详细信息<\/a>>\n","truncated":false}}
%---
%[output:04b8accc]
%   data: {"dataType":"textualVariable","outputData":{"name":"x","value":"-2.1500e-12"}}
%---
%[output:471c8b91]
%   data: {"dataType":"textualVariable","outputData":{"name":"fval","value":"-3.4106e-13"}}
%---
%[output:5fe99f88]
%   data: {"dataType":"text","outputData":{"text":"\n<a href = \"matlab: helpview('optim','eqn_solved','CSHelpWindow');\">方程已解<\/a>。\n\nfsolve 已完成，因为按照<a href = \"matlab: helpview('optim','fcn_tolerance_fsolve','CSHelpWindow');\">函数容差<\/a>的值衡量，\n函数值向量接近于零，并且按照梯度的值衡量，\n<a href = \"matlab: helpview('optim','appears_regular','CSHelpWindow');\">问题似乎为正则问题<\/a>。\n\n<<a href = \"matlab: createExitMsg({'optim:fsolve:Exit1basic','fsolve'},{'optim:fsolve:Exit1detailed','1.806021e-08','1.000000e-06','2.257309e-15','1.000000e-03'},true,true);;\">停止条件详细信息<\/a>>\n","truncated":false}}
%---
%[output:1a5a592c]
%   data: {"dataType":"textualVariable","outputData":{"name":"x","value":"-1.7407e-07"}}
%---
%[output:20b16d5a]
%   data: {"dataType":"textualVariable","outputData":{"name":"fval","value":"-4.7511e-08"}}
%---
%[output:5b7d4a6e]
%   data: {"dataType":"text","outputData":{"text":"\n<a href = \"matlab: helpview('optim','eqn_solved','CSHelpWindow');\">方程已解<\/a>。\n\nfsolve 已完成，因为按照<a href = \"matlab: helpview('optim','fcn_tolerance_fsolve','CSHelpWindow');\">函数容差<\/a>的值衡量，\n函数值向量接近于零，并且按照梯度的值衡量，\n<a href = \"matlab: helpview('optim','appears_regular','CSHelpWindow');\">问题似乎为正则问题<\/a>。\n\n<<a href = \"matlab: createExitMsg({'optim:fsolve:Exit1basic','fsolve'},{'optim:fsolve:Exit1detailed','1.418990e-10','1.000000e-06','8.198422e-21','1.000000e-03'},true,true);;\">停止条件详细信息<\/a>>\n","truncated":false}}
%---
%[output:1a112624]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"x","rows":1,"type":"double","value":[["0.1592","-1.4516"]]}}
%---
