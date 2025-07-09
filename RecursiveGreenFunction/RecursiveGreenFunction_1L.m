function G_1L = RecursiveGreenFunction_1L(HCC, V, layer, omega, Sigma_L, Sigma_R, eta, disorder_type, Gamma, Norb)
% 高效递归格林函数计算，支持多种无序类型
% HCC: 单层或多层哈密顿量
% V: 层间跃迁矩阵
% layer: 层数（若HCC为多层，自动识别）
% omega: 能量
% Sigma_L, Sigma_R: 左/右自能
% varargin{1}: eta（虚部，默认1e-5）
% varargin{2}: 'onsite_disorder'/'UnitCell_disorder'/'spin_disorder'
% varargin{3}: Gamma（无序强度）
% varargin{4}: N_orbital（单胞内轨道数，UnitCell无序用）
% 2025年6月24日

% %% 参数声明
% arguments (Input)
%     HCC (:,:) {mustBeNumeric}
%     V (:,:) {mustBeNumeric}
%     layer (1,1) {mustBeScalarOrEmpty,mustBeInteger}
%     omega (1,1) {mustBeNumeric}
%     Sigma_L (:,:) {mustBeNumeric}
%     Sigma_R (:,:) {mustBeNumeric}
%     eta (1,1) {mustBeNumeric}
%     disorder_type
%     Gamma (1,1) {mustBeNumeric}
%     Norb (1,1) {mustBeInteger}
% end
% 
% arguments (Output)
%     G_1L (:,:) {mustBeNumeric}
% end


%%
% RecursiveGreenFunction_1L
% Calculate the Retarded Greenfunction G_1L recursively HCC is the central
% area layer hamiltonian V is the central area hopping hamiltonian omega
% is the energy Sigma_L/R is the self energy of left/right lead
% varargin{1} is the amplitude of infinitisimal part.
% varargin{2} is the disorder part. If onsite, add random disorder to onsite energy
% varargin{3} is the parameter used to control the disorder
% Parameter setup 参数设定
% 若输入参数含有自设定无穷小虚部，则采用所设eta，否则默认10^(-5)
% varargin{1} = eta;
% varargin{2} = 'onsite_disorder';
%          or ='UnitCell_disorder'
%          or ="spin_disorder"
% varargin{3} = Gamma; % Disorder
% varargin{4} = N_orbital % 单胞内轨道数
if isempty(eta)
    eta = 1e-5;
end
% 若不同层的哈密顿量是全同的，那么HCC应当是(W,W)型矩阵。
% 当不同层哈密顿量不是全同的，则HCC应当形如(W, W*Layer)式的矩阵。
[Width, Length] = size(HCC);
% 判定层的全同性。若不是全同的(layerflag =0)，那么每一层都应当单独迭代
layerflag = (Width==Length);
Iden = eye(Width);
M_ii = complex(zeros(Width));
G_1L = complex(zeros(Width));
% 循环体计算
% 当不同层的哈密顿量是全同的，那么HCC应当是(W,W)型矩阵。 跃迁项V也是(W,W)型矩阵。此时迭代层数需要手动输入定义为layer
% 判定层是否全同
if layerflag
    % 全同层情形
    % 判定是否加入单格点在位无序
    
    if strcmp(disorder_type,"Onsite_disorder")
        H_impurity = diag((rand(1,Width)-0.5)*Gamma);
        M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - Sigma_L) \ Iden;
        G_1L = M_ii;
        for ii = 2:layer
            H_impurity = diag((rand(1,Width)-0.5)*Gamma);
            M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - V'*M_ii*V) \ Iden;
            G_1L = G_1L * V * M_ii;
        end

        % 判断是否加入单胞无序
    elseif strcmp(disorder_type,"UnitCell_disorder")
        H_impurity = diag(repelem((rand(1, Width/Norb)-0.5)*Gamma, Norb));
        M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - Sigma_L) \ Iden;
        G_1L = M_ii;

        for ii = 2:layer
            H_impurity = diag(repelem((rand(1, Width/Norb)-0.5)*Gamma, Norb));
            M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - V'*M_ii*V) \ Iden;
            G_1L = G_1L * V * M_ii;
        end
    elseif strcmp(disorder_type,"spin_disorder")
        % 判定是否加入自旋无序
        sigma_x = [0,1;1,0];
        sigma_y = [0,-1i;1i,0];
        H_impurity = (diag(repelem((rand(1, Width/2)-0.5), 2)) +...
                      kron(diag(rand(1,Width/2)-0.5), sigma_x) +...
                      kron(diag(rand(1,Width/2)-0.5), sigma_y))*Gamma;
        M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - Sigma_L) \ Iden;
        G_1L = M_ii;
        for ii = 2:layer
            H_impurity = (diag(repelem((rand(1, Width/2)-0.5), 2)) +...
                          kron(diag(rand(1,Width/2)-0.5), sigma_x) +...
                          kron(diag(rand(1,Width/2)-0.5), sigma_y))*Gamma;
            M_ii = ((omega + 1i * eta)*Iden - HCC - H_impurity - V'*M_ii*V) \ Iden;
            G_1L = G_1L * V * M_ii;
        end
    elseif strcmp(disorder_type,"Clean")
        % 不加入无序
        M_ii = ((omega + 1i * eta)*Iden - HCC - Sigma_L) \ Iden;
        G_1L = M_ii;
        for ii = 2:layer
            M_ii = ((omega + 1i * eta)*Iden - HCC - V'*M_ii*V) \ Iden;
            G_1L = G_1L * V * M_ii;
        end
    end
    % 右自能修正
    M_ii = (M_ii\Iden - Sigma_R) \ Iden;
    G_1L = G_1L + G_1L * Sigma_R * M_ii;
else
    % 当不同层的哈密顿量不全同的，那么HCC应当是(W,W*Layer)型矩阵
    % 跃迁项V则是(W,W*(Layer-1))型矩阵。此时

    M_ii = ((omega + 1i * eta)*Iden - HCC(:, 1:Width) - Sigma_L) \ Iden;
    G_1L = M_ii;
    for ii = 2:Length/Width
        H_ii = HCC(:, ((ii-1)*Width+1): ii*Width);
        V_ii = V(:,((ii-2)*Width+1): (ii-1)*Width);
        M_ii = ((omega + 1i * eta)*Iden - H_ii - V_ii'*M_ii*V_ii) \ Iden;
        G_1L = G_1L * V_ii * M_ii;
    end
    M_ii = (M_ii\Iden - Sigma_R) \ Iden;
    G_1L = G_1L + G_1L * Sigma_R * M_ii;
end
end