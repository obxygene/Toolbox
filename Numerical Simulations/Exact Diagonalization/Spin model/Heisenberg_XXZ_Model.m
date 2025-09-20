% 主脚本：构造并显示Heisenberg XXZ模型的哈密顿量
clear; % 清除工作区变量
L = 10;           % 链长
J = 1.0;          % Heisenberg相互作用强度
delta = 0;        % Ising各向异性参数
sz = 2;           % 总磁量子数Sz
k = 1;            % 动量子空间编号

% 计算哈密顿量的稀疏矩阵三元组
[row, col, data] = get_hamiltonian_sparse(L, J, delta, sz, k);
sizeH = max([row, col]) + 1; % 矩阵维度
H = sparse(row+1, col+1, data, sizeH, sizeH); % MATLAB索引从1开始
disp(full(H)) % 显示稠密矩阵

%----------------- 主函数及子函数 -----------------%

function [hamiltonian_rows, hamiltonian_cols, hamiltonian_data] = get_hamiltonian_sparse(L, J, delta, sz, k)
    % 构造Heisenberg XXZ模型的哈密顿量，动量守恒
    % 输入:
    %   L     - 链长
    %   J     - Heisenberg相互作用强度
    %   delta - Ising各向异性参数
    %   sz    - 总磁量子数Sz
    %   k     - 动量子空间编号
    % 输出:
    %   hamiltonian_rows, hamiltonian_cols, hamiltonian_data - 稀疏矩阵三元组

    % 检查Sz取值是否合法
    assert(sz <= (floor(L/2) + mod(L,2)) && sz >= -floor(L/2));

    % 生成动量本征态基矢及其归一化因子
    basis_states = [];
    norms = [];
    state = first_state(L, sz);
    end_state = last_state(L, sz);
    while state <= end_state
        [representative, translation] = get_representative(L, state);
        if state == representative
            amplitude = 0.0;
            for n_translation_sites = 0:L-1
                new_state = translate(L, state, n_translation_sites);
                if new_state == state
                    amplitude = amplitude + exp(1i*2*pi*k/L*n_translation_sites); % 动量因子
                end
            end
            norm = sqrt(abs(amplitude));
            if norm > 1e-12
                basis_states(end+1) = state;
                norms(end+1) = norm;
            end
        end
        state = next_state(state);
    end

    % 构造周期性链的键
    heisenberg_bonds = [ (0:L-1)' mod((0:L-1)'+1, L) ];
    n_states = length(basis_states); % 基矢数量
    n_bonds = size(heisenberg_bonds,1); % 键数量
    max_nnz = n_states * (1 + 2*n_bonds); % 预分配最大非零元数
    hamiltonian_rows = zeros(1, max_nnz);
    hamiltonian_cols = zeros(1, max_nnz);
    hamiltonian_data = zeros(1, max_nnz);
    nnz = 0;

    % 遍历所有基矢，填充哈密顿量的非零元
    for state_index = 1:n_states
        state = basis_states(state_index);
        diagonal = 0;
        % 计算对角项（Ising项）
        for b = 1:n_bonds
            i = heisenberg_bonds(b,1);
            j = heisenberg_bonds(b,2);
            if get_site_value(state, i) == get_site_value(state, j)
                diagonal = diagonal + J/4;
            else
                diagonal = diagonal - J/4;
            end
        end
        nnz = nnz + 1;
        hamiltonian_rows(nnz) = state_index-1;
        hamiltonian_cols(nnz) = state_index-1;
        hamiltonian_data(nnz) = diagonal;
        % 计算非对角项（自旋交换）
        for b = 1:n_bonds
            i = heisenberg_bonds(b,1);
            j = heisenberg_bonds(b,2);
            if get_site_value(state, i) ~= get_site_value(state, j)
                flipmask = bitshift(1, i) + bitshift(1, j);
                new_state = bitxor(state, flipmask);
                [representative, translation] = get_representative(L, new_state);
                new_state_index = find(basis_states == representative, 1) - 1;
                coeff = ((J+delta)/2)*(norms(new_state_index+1)/norms(state_index))*exp(1i*2*pi*k/L*translation)/2;
                nnz = nnz + 1;
                hamiltonian_rows(nnz) = new_state_index;
                hamiltonian_cols(nnz) = state_index-1;
                hamiltonian_data(nnz) = coeff;
                nnz = nnz + 1;
                hamiltonian_rows(nnz) = state_index-1;
                hamiltonian_cols(nnz) = new_state_index;
                hamiltonian_data(nnz) = conj(coeff);
            end
        end
    end
    hamiltonian_rows = hamiltonian_rows(1:nnz);
    hamiltonian_cols = hamiltonian_cols(1:nnz);
    hamiltonian_data = hamiltonian_data(1:nnz);
end

% 获取某一位的自旋（0/1）
function val = get_site_value(state, site)
    val = bitget(state, site+1);
end

% 设置某一位的自旋值
function state = set_site_value(state, site, value)
    mask = bitshift(1, site);
    if value
        state = bitor(state, mask);
    else
        state = bitand(state, bitcmp(mask, 'int64'));
    end
end

% 返回第一个满足Sz的基态
function state = first_state(L, sz)
    n_upspins = floor(L/2) + sz;
    state = bitshift(1, n_upspins) - 1;
end

% 生成下一个具有相同上自旋数的状态（二进制排列）
function state_out = next_state(state)
    t = bitor(state, state - 1, 'int64') + 1;
    numer = bitand(t, -t, 'int64');
    denom = bitand(state, -state, 'int64');
    state_out = bitor(t, bitshift(floorDiv(numer, denom), -1, 'int64') - 1, 'int64');
end

% 返回最后一个满足Sz的基态
function state = last_state(L, sz)
    n_upspins = floor(L/2) + sz;
    state = bitshift((bitshift(1, n_upspins) - 1), L - n_upspins);
end

% 对状态进行平移
function new_state = translate(L, state, n_translation_sites)
    new_state = int64(0);
    for site = 0:L-1
        site_value = get_site_value(state, site);
        new_state = set_site_value(new_state, mod(site + n_translation_sites, L), site_value);
    end
end

% 找到状态的最小表示及其平移量
function [representative, translation] = get_representative(L, state)
    representative = state;
    translation = 0;
    for n_translation_sites = 0:L-1
        new_state = translate(L, state, n_translation_sites);
        if new_state < representative
            representative = new_state;
            translation = n_translation_sites;
        end
    end
end %[output:2f80d563]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:2f80d563]
%   data: {"dataType":"text","outputData":{"text":"   1.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.4045 - 0.2939i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i\n   0.5000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.1545 - 0.4755i   0.4045 + 0.2939i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.4045 + 0.2939i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.4045 + 0.2939i   0.0000 + 0.0000i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.4045 + 0.2939i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i  -0.1545 + 0.4755i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i\n   0.4045 + 0.2939i   0.1545 + 0.4755i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i   0.5000 + 0.0000i   0.1545 + 0.4755i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.4045 - 0.2939i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.1545 - 0.4755i   0.1545 - 0.4755i  -0.5000 + 0.0000i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.1545 - 0.4755i   0.0000 + 0.0000i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.4045 - 0.2939i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i  -0.5000 + 0.0000i   0.5000 + 0.0000i  -0.1545 - 0.4755i   0.4045 + 0.2939i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.4045 - 0.2939i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.5000 + 0.0000i  -0.5000 + 0.0000i   0.5000 + 0.0000i   0.0955 + 0.2939i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.4045 - 0.2939i   0.5000 + 0.0000i   0.0000 + 0.0000i   0.1545 + 0.4755i  -0.1545 + 0.4755i   0.5000 + 0.0000i  -0.5000 + 0.0000i  -0.1545 + 0.4755i\n   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.4045 - 0.2939i   0.0955 - 0.2939i  -0.1545 - 0.4755i  -0.8090 + 0.0000i\n\n","truncated":false}}
%---
