function H = Graphene_Armchair_Real_Sparse_V2(Length, Width, t, options)
%% Graphene nanoribbon with armchair edges
%     1     2               A   B           w=1
%   3         4         B           A
%     5     6               A   B           w=2
%       ...                                 ...
%  4W-3    4W-2             A   B           w=Width
% 4W-1        4W        B               A
%   4W+1   4W+2            A   B           w=Width+1
% --------------------------------------
% 4W+3        4W+4      B           A       w=Width+1
%
% Usage:
%   H = Graphene_Armchair_Real_Sparse_1(10, 5, -2.7)
%   H = Graphene_Armchair_Real_Sparse_1(10, 5, -2.7, Boundary="PBC")
%   H = Graphene_Armchair_Real_Sparse_1(10, 5, -2.7, Boundary="xPBC")
%   H = Graphene_Armchair_Real_Sparse_1(10, 5, -2.7, Bdirection="Bx", Bfield=0.1)
%   H = Graphene_Armchair_Real_Sparse_1(10, 5, -2.7, Bdirection="By", Bfield=0.2, Boundary="PBC")

% 参数验证块
arguments
    Length (1,1) {mustBePositive, mustBeInteger}
    Width (1,1) {mustBePositive, mustBeInteger}
    t (1,1) {mustBeNumeric}
    options.Boundary {mustBeText, mustBeMember(options.Boundary, ["PBC", "OBC", "xPBC", "yPBC"])} = "OBC"
    options.Bdirection {mustBeText, mustBeMember(options.Bdirection, ["Bx", "By", "" ] )} = ""
    options.Bfield (1,1) {mustBeNumeric} = 0
end

% 解析边界条件
switch options.Boundary
    case "PBC"
        boundary_x = "PBC";
        boundary_y = "PBC";
    case "OBC"
        boundary_x = "OBC";
        boundary_y = "OBC";
    case "xPBC"
        boundary_x = "PBC";
        boundary_y = "OBC";
    case "yPBC"
        boundary_x = "OBC";
        boundary_y = "PBC";
end

% 基本哈密顿量矩阵
H_intra_00 = [0, t, t, 0;
              t, 0, 0, t;
              t, 0, 0, 0;
              0, t, 0, 0;];

H_intra_0y = [0, 0, 0, 0;
              0, 0, 0, 0;
              t, 0, 0, 0;
              0, t, 0, 0;];
H_intra_y0 = H_intra_0y';

% 构建跃迁矩阵
[Hopping_x, Hopping_y, total_atoms, V, y_matrix_size] = buildHoppingMatrices(...
    Length, Width, boundary_x, boundary_y, t);

% 根据磁场条件构建哈密顿量
if ~isempty(options.Bdirection) && options.Bfield ~= 0
    H = buildHamiltonianWithField(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
        Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, ...
        options.Bdirection, options.Bfield);
else
    H = buildHamiltonianNoField(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
        Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size);
end
end

%% 辅助函数：构建跃迁矩阵
function [Hopping_x, Hopping_y, total_atoms, V, y_matrix_size] = buildHoppingMatrices(...
    Length, Width, boundary_x, boundary_y, t)

    % x方向跃迁矩阵（尺寸满足 spdiags 约束）
    Hopping_x = spdiags(ones(Length,1), 1, Length, Length);
    if boundary_x == "PBC"
        Hopping_x(Length, 1) = 1;
    end

    % y方向跃迁矩阵和原子总数（统一为 (Width+1)x(Width+1)）
    y_matrix_size = Width + 1;
    if boundary_y == "PBC"
        % 注意：spdiags(B,d,m,n) 要求 B 为 length min(m,n) 的列向量
        Hopping_y = spdiags(ones(y_matrix_size,1), 1, y_matrix_size, y_matrix_size);
        Hopping_y(y_matrix_size, 1) = 1;  % 回绕
        total_atoms = 4 * y_matrix_size;  % 4*(Width+1)
    else
        Hopping_y = spdiags(ones(y_matrix_size,1), 1, y_matrix_size, y_matrix_size);
        total_atoms = 4 * Width + 2;      % OBC 需裁剪
    end

    % 层间跃迁矩阵
    V_inter = zeros(4,4);
    V_inter(4,3) = t;

    % V 的组装
    V = kron(eye(y_matrix_size), V_inter);
    if boundary_y == "OBC"
        V = sparse(V(1:total_atoms, 1:total_atoms));
    else
        V = sparse(V);
    end
end

%% 辅助函数：无磁场哈密顿量
function H = buildHamiltonianNoField(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
    Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size)

    % y 方向层内与层间（保证三项尺寸一致）
    H = kron(eye(y_matrix_size), H_intra_00) ...
      + kron(Hopping_y,  H_intra_0y) ...
      + kron(Hopping_y', H_intra_y0);

    % OBC 时裁剪掉 armchair 最后一层多余原子
    if boundary_y == "OBC"
        H = sparse(H(1:total_atoms, 1:total_atoms));
    else
        H = sparse(H);
    end

    % x 方向层之间的耦合
    H = kron(speye(Length), H) + kron(Hopping_x, V) + kron(Hopping_x', V');
end

%% 辅助函数：有磁场哈密顿量总控
function H = buildHamiltonianWithField(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
    Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, Bdirection, Bfield)

    % 使用下三角确保后续 temp + temp' 保持厄米
    H_intra_00_tril = tril(H_intra_00);

    if Bdirection == "Bx"
        H = buildHamiltonianWithBx(Length, Width, H_intra_00_tril, H_intra_0y, H_intra_y0, ...
            Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, Bfield);

    elseif Bdirection == "By"
        H = buildHamiltonianWithBy(Length, Width, H_intra_00_tril, H_intra_0y, H_intra_y0, ...
            Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, Bfield);
    else
        % 兜底：方向无效则退化为无磁场
        H = buildHamiltonianNoField(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
            Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size);
    end
end

%% Bx 磁场处理
function H = buildHamiltonianWithBx(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
    Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, Bfield)

    % 基块（使用下三角，已由上层传入）
    H_base  = kron(eye(y_matrix_size), H_intra_00);
    temp_hop = kron(Hopping_y, H_intra_0y);

    % OBC 时裁剪掉 armchair 最后一层多余原子
    if boundary_y == "OBC"
        H_base   = sparse(H_base(1:total_atoms, 1:total_atoms));
        temp_hop = sparse(temp_hop(1:total_atoms, 1:total_atoms));
    else
        H_base   = sparse(H_base);
        temp_hop = sparse(temp_hop);
    end

    % 相位矩阵（修正尺寸）
    D_ones  = speye(Length);
    D_field = spdiags(exp(1i*(0:Length-1)*Bfield)', 0, Length, Length);

    temp = kron(D_ones,  H_base) ...
         + kron(D_field, temp_hop);

    % 组合得到整条带的哈密顿量（厄米化 + x 向耦合）
    H = temp + temp' + kron(Hopping_x, V) + kron(Hopping_x', V');
end

%% By 磁场处理
function H = buildHamiltonianWithBy(Length, Width, H_intra_00, H_intra_0y, H_intra_y0, ...
    Hopping_x, Hopping_y, V, total_atoms, boundary_y, y_matrix_size, Bfield)

    phase_vector = exp(1i * (1:y_matrix_size) * Bfield);
    phase_matrix = diag(phase_vector);

    temp    = kron(phase_matrix, H_intra_00);
    H_layer = temp + temp' + kron(Hopping_y, H_intra_0y) + kron(Hopping_y', H_intra_y0);

    % OBC 时裁剪
    if boundary_y == "OBC"
        H_layer = sparse(H_layer(1:total_atoms, 1:total_atoms));
    else
        H_layer = sparse(H_layer);
    end

    % x 方向层之间的耦合
    H = kron(speye(Length), H_layer) + kron(Hopping_x, V) + kron(Hopping_x', V');
end
