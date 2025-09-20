clear;
L = 6; J = 1.0; hs = 7; sz = 1;        % 总 Sz = 0 子空间
H = full(get_hamiltonian_sparse(L, J, hs, sz));
disp(H) %[output:867e0c2a]
% eigs_small = eig(full(H));               % 小尺寸可直接求全谱
%[text] $H=J\\sum\_{\\langle{i,j}\\rangle}\\mathbf{S}\_{i}\\cdot\\mathbf{S}\_{j}+h\_{S}\\sum\_{i}(-1)^{i}S^{z}\_i$
function H = get_hamiltonian_sparse(L, J, hs, sz)
% Function to create the Hamiltonian of the spin 1/2 Heisenberg model
% with staggered magnetic field, Sz conservation implemented
%
% Args:
%   L  (int): length of chain
%   J  (float): coupling constant for Heisenberg term
%   hs (float): coupling constant for staggered field
%   sz (int): total Sz
%
% Returns:
%   hamiltonian_rows (array): row index of non-zero elements
%   hamiltonian_cols (array): column index of non-zero elements
%   hamiltonian_data (array): value of non-zero elements



    function state = first_state(L, sz)
        n_upspins = floorDiv(L,2) + sz;
        state = bitshift(1, n_upspins) - 1;
    end

    function state_out = next_state(state)
        t = bitor(state, state - 1, 'int64') + 1;
        numer = bitand(t, -t, 'int64');
        denom = bitand(state, -state, 'int64');
        state_out = bitor(t, bitshift(floorDiv(numer, denom), -1, 'int64') - 1, 'int64');
    end

    function state = last_state(L, sz)
        n_upspins = floorDiv(L,2) + sz;
        state = bitshift((bitshift(1, n_upspins) - 1), L - n_upspins);
    end

% Check if sz is valid
assert(sz <= (floorDiv(L,2) + mod(L,2)) && sz >= -floorDiv(L,2));

% Create list of states with fixed sz
basis_states = [];
state = first_state(L, sz);
end_state = last_state(L, sz);
while state <= end_state
    basis_states(end+1) = state; %#ok<AGROW>
    state = next_state(state);
end

% Define chain lattice
heisenberg_bonds = [ (0:L-1)' mod((0:L-1)'+1, L) ];


% 预分配稀疏矩阵空间
n_states = length(basis_states);
n_bonds = size(heisenberg_bonds,1);
max_nnz = n_states * (1 + n_bonds); % 每个状态最多1个对角+每个bond一个非对角
hamiltonian_rows = zeros(1, max_nnz);
hamiltonian_cols = zeros(1, max_nnz);
hamiltonian_data = zeros(1, max_nnz);
nnz = 0;

% Run through all spin configurations
for state_index = 1:n_states
    state = basis_states(state_index);

    % Apply diagonal Ising bonds
    diagonal = 0;
    for b = 1:n_bonds
        i = heisenberg_bonds(b,1);
        j = heisenberg_bonds(b,2);
        if bitget(state, i+1) == bitget(state, j+1)
            diagonal = diagonal + J/4;
        else
            diagonal = diagonal - J/4;
        end
        % disp(diagonal)
    end

    % Apply diagonal staggered Sz field
    for site = 0:2:L-1
        diagonal = diagonal + hs * (2*bitget(state, site+1) - 1);
        diagonal = diagonal - hs * (2*bitget(state, site+2) - 1);
    end

    nnz = nnz + 1;
    hamiltonian_rows(nnz) = state_index;
    hamiltonian_cols(nnz) = state_index;
    hamiltonian_data(nnz) = diagonal;

    % Apply exchange interaction
    for b = 1:n_bonds
        i = heisenberg_bonds(b,1);
        j = heisenberg_bonds(b,2);
        if bitget(state, i+1) ~= bitget(state, j+1)
            flipmask = bitor(bitshift(1, i), bitshift(1, j));
            new_state = bitxor(state, flipmask);
            new_state_index = find(basis_states == new_state, 1);
            nnz = nnz + 1;
            hamiltonian_rows(nnz) = state_index;
            hamiltonian_cols(nnz) = new_state_index;
            hamiltonian_data(nnz) = J/2;
        end
    end
end
% 截断到实际长度
hamiltonian_rows = hamiltonian_rows(1:nnz);
hamiltonian_cols = hamiltonian_cols(1:nnz);
hamiltonian_data = hamiltonian_data(1:nnz);
H = sparse(hamiltonian_rows, hamiltonian_cols, hamiltonian_data);
end





%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:867e0c2a]
%   data: {"dataType":"text","outputData":{"text":"    0.5000    0.5000         0         0         0         0         0         0    0.5000         0         0         0         0         0         0\n    0.5000   27.5000    0.5000         0         0    0.5000         0         0         0         0         0    0.5000         0         0         0\n         0    0.5000   -0.5000    0.5000         0         0    0.5000         0         0         0         0         0         0    0.5000         0\n         0         0    0.5000   27.5000    0.5000         0         0    0.5000         0         0         0         0         0         0    0.5000\n         0         0         0    0.5000    0.5000         0         0         0    0.5000         0         0         0         0         0         0\n         0    0.5000         0         0         0    0.5000    0.5000         0         0         0         0         0         0         0         0\n         0         0    0.5000         0         0    0.5000  -28.5000    0.5000         0    0.5000         0         0         0         0         0\n         0         0         0    0.5000         0         0    0.5000   -0.5000    0.5000         0    0.5000         0         0         0         0\n    0.5000         0         0         0    0.5000         0         0    0.5000  -28.5000         0         0    0.5000         0         0         0\n         0         0         0         0         0         0    0.5000         0         0    0.5000    0.5000         0         0         0         0\n         0         0         0         0         0         0         0    0.5000         0    0.5000   27.5000    0.5000    0.5000         0         0\n         0    0.5000         0         0         0         0         0         0    0.5000         0    0.5000   -0.5000         0    0.5000         0\n         0         0         0         0         0         0         0         0         0         0    0.5000         0    0.5000    0.5000         0\n         0         0    0.5000         0         0         0         0         0         0         0         0    0.5000    0.5000  -28.5000    0.5000\n         0         0         0    0.5000         0         0         0         0         0         0         0         0         0    0.5000    0.5000\n\n","truncated":false}}
%---
