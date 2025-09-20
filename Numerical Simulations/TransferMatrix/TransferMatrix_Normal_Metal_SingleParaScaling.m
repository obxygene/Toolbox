clear;
rng(1)
% Define a marker database
markers = {'o', 's', 'd', '*', '^', 'p', 'h'}; % Circle, square, diamond, star, triangle, pentagon, hexagon
colors = {'r', 'g', 'b', 'm', 'c', 	"#EDB120", 'k',"#A2142F"}; % Red, green, blue, magenta, cyan, yellow, black
tic;

Central.N_UnitCell = 1;
Central.N_row = 2; %[control:editfield:7718]{"position":[17,18]}
Central.N_col = 3000; %[control:editfield:6391]{"position":[17,21]}
BoundaryCond = "OBC"; %[control:dropdown:4745]{"position":[16,21]}
Central.BoundaryCond = BoundaryCond;
epsilon = 0 * ones(1,Central.N_UnitCell);
Central.mu = 0;
% 中心区其他量
% 总格点数
N_tot_site = Central.N_row * Central.N_col * Central.N_UnitCell;
N_Cells = Central.N_row * Central.N_col;
% 每列格点数
N_site_pcol = Central.N_row * Central.N_UnitCell;
% 跃迁
t = -1;
% Fermi Energy
FermiEnergy = 0;
% Disorder samplings
N_Disorder = 11; %[control:editfield:3d70]{"position":[14,16]}
Disorder = linspace(0,10,N_Disorder);
% Quasi 1D Width
Width = 6:2:30;


LocLength = zeros(length(Width), length(Disorder));
% Loop for disorder and delta simutaneously
for jj = 1:length(Width)
    Central.N_row = Width(jj);
    [H_CC, V_CC] = Hamiltonian_Metal(epsilon, t, Central.N_row);
    r = rank(V_CC);
    R = eye(size(H_CC,1));
    R_Spectrum = zeros(2*r, Central.N_col);
    for kk = 1:length(Disorder)
        for ii = 1:Central.N_col
            if ii > 1
                T = Transfer_Matrix(H_CC, V_CC, FermiEnergy, Disorder(kk))*Q;
            else
                T = Transfer_Matrix(H_CC, V_CC, FermiEnergy, Disorder(kk));
            end
            [Q, R] = qr(T);
            [R_diag, ind] = sort(diag(R), 'descend','ComparisonMethod','abs');
            Q = Q(:, ind);
            R_Spectrum(:, ii) = abs(R_diag);
        end
        lnR = log(R_Spectrum);        
        Lyapunov_Spectrum = sum(lnR,2)/Central.N_col;        
    % The minimal Lyapunov exponent corresponds to the largest localization
    % length
        LocLength(jj, kk) = 1./min(abs(Lyapunov_Spectrum));
    end
end
%%
figure
plot(Disorder, LocLength./Width','o-');
xlabel('Disorder')
ylabel('\xi/M')
legend(arrayfun(@(x) ['M =' num2str(x)], Width, 'UniformOutput',false))
yscale log
% xscale log
hold off

Lambda = LocLength ./ Width';
%[text] 利用拟合得到$\\xi\_\\infty$
figure; xscale log; yscale log; hold on

eqn = 'log(1+k*xi/x)/k'; startpoint = [1,1]; Lower = [0,0];

xi_infty = zeros(1,length(Disorder));
xi_over_Width = zeros(length(Width), length(Disorder));

for ii = 1:length(Disorder)
    ind = ~isnan(Lambda(:, ii));
    f = fit(Width(ind(4:end))', Lambda(ind(4:end), ii), eqn, "startpoint",startpoint,"lower",Lower);
    xi_infty(ii) = f.xi;
    xi_over_Width(:, ii) = xi_infty(ii)./Width;    
    PlotFlag = "Fit"; %[control:dropdown:597e]{"position":[16,21]}
    if PlotFlag == "Fit"
        xscale log; yscale log
        p1 = plot(Width, Lambda(:, ii),'o','DisplayName',['W=',num2str(Disorder(ii))]);
        p2 = plot(Width, f(Width),'-','LineWidth',2,'DisplayName',['\xi=',num2str(f.xi)]);
        xlabel('$M$','Interpreter','latex')
        ylabel('$\Lambda$','Interpreter','latex')
        % pause 
    else
        % Scaling plot
        xi_infty(ii) = f.xi; % Fit parameter of xi_infty for a given width
        p = plot(xi_over_Width(:, ii), Lambda(:,ii),'.','DisplayName',['\Gamma=',num2str(Disorder(ii))]);
        p.LineStyle = '-';
        p.Marker = markers{mod(ii,length(markers))+1};
        p.MarkerEdgeColor=colors{mod(ii,length(colors))+1};
        p.MarkerSize = 6;
    end
end
legend

%%
% % [xi_d_Width, ind] = sort(xi_d_Width(:));
% % Lambda_sort = Lambda(ind);
% % ind = (~isnan(Lambda(:))) & (xi_d_Width > 0);
% % figure;
% ind = ~isnan(Lambda(:)) & (xi_over_Width(:)<0.1);
% Width_Mat = repmat(Width, length(Disorder));
% startpoint = [1,1];
% Lower = [0,0];
% eqn = 'xi/x-a*(xi/x)^2';
% f2 = fit(Width_Mat(ind), Lambda(ind), eqn, "startpoint",startpoint,"lower",Lower);
% % p = plot([1e-5;xi_d_Width(ind)], f2([1e-5;xi_d_Width(ind)]));
% legend('AutoUpdate','off')
% plot((f.xi)./logspace(-15,3,101), f(logspace(-15,3,101)),'Color','black','linewidth',2);
% 
% xlabel('$\log(\xi_\infty/M)$','Interpreter','latex')
% ylabel('$\log\Lambda$','Interpreter','latex')
% 
% box on














%%
%[text] ### 辅助函数
toc; %[output:5caee14d]
function T = Transfer_Matrix(H_CC, V_CC, FermiEnergy, Disorder)
eta = 0;
[V, Xi, W] = svd(V_CC);
r = rank(Xi);
V = V(:, 1:r);
W = W(:, 1:r);
Xi = Xi(1:r, 1:r);
H_imp = Disorder * diag(rand(1,size(H_CC,1))-0.5);
G = inv((FermiEnergy + 1i * eta) * eye(size(H_CC,1)) - H_CC - H_imp);
G_vv = (V' * G * V) * Xi;
G_ww = (W' * G * W) * Xi;
G_vw = (W' * G * V) * Xi;
G_wv = (V' * G * W) * Xi;

theta = angle(det(G_vw));


T = [inv(G_vw), -inv(G_vw) * G_ww; 
    G_vv * inv(G_vw), G_wv - G_vv * inv(G_vw) * G_ww];
T = exp(1i * theta / r) * T;
end


function [H00, H01] = Hamiltonian_Metal(epsilon, t, Width)
    H00 = epsilon * eye(Width) + t * (diag(ones(1, Width-1), +1) + diag(ones(1, Width-1), -1));
    H01 = t * eye(Width);
end



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[control:editfield:7718]
%   data: {"defaultValue":0,"label":"编辑字段","run":"Section","valueType":"Double"}
%---
%[control:editfield:6391]
%   data: {"defaultValue":0,"label":"编辑字段","run":"Section","valueType":"Double"}
%---
%[control:dropdown:4745]
%   data: {"defaultValue":"\"OBC\"","itemLabels":["OBC","PBC"],"items":["\"OBC\"","\"PBC\""],"label":"下拉列表","run":"Section"}
%---
%[control:editfield:3d70]
%   data: {"defaultValue":0,"label":"编辑字段","run":"Section","valueType":"Double"}
%---
%[control:dropdown:597e]
%   data: {"defaultValue":"\"Fit\"","itemLabels":["拟合结果","单参数标度结果"],"items":["\"Fit\"","\"SingleParam\""],"label":"PlotFlat","run":"AllSections"}
%---
%[output:5caee14d]
%   data: {"dataType":"text","outputData":{"text":"历时 86.248624 秒。\n","truncated":false}}
%---
