function H = Graphene_ZigZag_Real_Sparse(Length, Width, t, varargin)
%% labeling
% 1 -- 2                    B -- A      w = 1
%      |                         |
% 4 -- 3                    A -- B
% |                         |
% 5 -- 6                    B -- A      w = 2
%      |
% 8 -- 7
% |
%  ....
% 2W-3 -- 2W-2              B -- A      w = Width
%           |                    |
%  2W  -- 2W-1              A -- B
%
% 单胞内跃迁
% varargin{1} [string] 选择(0, Bx, 0)或(By, 0, 0)
%               "Bx" or "By"
% varargin{2} [scalar] 给出磁场
% 单胞内跃迁


H_intra_00 = [0, t, 0, 0;
    t, 0, t, 0;
    0, t, 0, t;
    0, 0, t, 0;];
% 单胞间跃迁
H_intra_0y = [0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    t, 0, 0, 0;];
H_intra_y0 = H_intra_0y';

Hopping_x = spdiags(ones(Length,1),1,Length,Length);
Hopping_y = diag(ones(1,Width-1),1);

% 构建层间跃迁
V_inter = [0, 0, 0, 0;
    t, 0, 0, 0;
    0, 0, 0, t;
    0, 0, 0, 0];
V = kron(eye(Width), V_inter);
V = sparse(V);


if nargin > 3
    H_intra_00 = tril(H_intra_00);
    if varargin{1} == "(0,Bx)"

        H = sparse(kron(eye(Width), H_intra_00));
        temp = sparse(kron(Hopping_y, H_intra_0y));

        temp = kron(spdiags(ones(Length,1),0,Length,Length),H) + ...
            kron(spdiags(exp(1i*(0:Length)*varargin{2})',0,Length,Length), temp);
        H = temp + temp' + kron(Hopping_x,V) + kron(Hopping_x',V');

    elseif varargin{1} == "(-By,0)"

        temp = kron(diag(exp(1i*(1:Width)*varargin{2})), H_intra_00);
        H = temp + temp' + kron(Hopping_y, H_intra_0y)...
            + kron(Hopping_y', H_intra_y0);

        H = sparse(H);
        H = kron(speye(Length), H) + kron(Hopping_x,V) + kron(Hopping_x',V');
    end
else
    % 构成单层的哈密顿量
    H = kron(eye(Width), H_intra_00) + kron(Hopping_y, H_intra_0y) + kron(Hopping_y', H_intra_y0);
    H = sparse(H);
    % 构成层间整体哈密顿量
    H = kron(speye(Length), H) + kron(Hopping_x, V) + kron(Hopping_x', V');
end


end