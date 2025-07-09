function H = Graphene_Armchair_Real_Sparse(Length, Width, t, varargin)
%% labeling
%     1     2               A   B           w=1
%   3         4         B           A
%     5     6               A   B           w=2
%       ...                                 ...
%  4W-3    4W-2             A   B           w=Width
% 4W-1        4W        B           A
%   4W+1   4W+2            A   B           w=Width+1
% --------------------------------------
% 4W+1        4W+2      B           A       w=Width+1
% 单胞内跃迁
% varargin{1} [string] 选择(0, Bx, 0)或(By, 0, 0)
%               "Bx" or "By"
% varargin{2} [scalar] 给出磁场

% 单胞内跃迁
H_intra_00 = [0, t, t, 0;
    t, 0, 0, t;
    t, 0, 0, 0;
    0, t, 0, 0;];

% 单胞间跃迁
H_intra_0y = [0, 0, 0, 0;
    0, 0, 0, 0;
    t, 0, 0, 0;
    0, t, 0, 0;];
H_intra_y0 = H_intra_0y';

Hopping_x = spdiags(ones(Length,1),1,Length,Length);
Hopping_y = diag(ones(1,Width),1);

% 构建层间跃迁
V_inter = zeros(4,4);
V_inter(4,3) = t;
V = kron(eye(Width+1), V_inter);
V = sparse(V(1:4*Width+2, 1:4*Width+2));


if nargin > 3
    H_intra_00 = tril(H_intra_00);
    if varargin{1} == "Bx"

        % H = kron(eye(Width+1), H_intra_00) + kron(Hopping_y, H_intra_0y);
        % % 由于armchair的石墨烯中，最后一层原子不存在，
        % % 构建出来的哈密顿量需要删去单胞构造多出来的一层原子
        % H = sparse(H(1:4*Width+2, 1:4*Width+2));
        % temp = kron(spdiags(exp(1i*(0:Length)*varargin{2})',0,Length,Length), H);
        % H = temp + temp' + kron(Hopping_x,V) + kron(Hopping_x',V');

        H = kron(eye(Width+1), H_intra_00);
        temp = kron(Hopping_y, H_intra_0y);
        
        % 由于armchair的石墨烯中，最后一层原子不存在，
        % 构建出来的哈密顿量需要删去单胞构造多出来的一层原子
        H = sparse(H(1:4*Width+2, 1:4*Width+2));
        temp = sparse(temp(1:4*Width+2, 1:4*Width+2));
        temp = kron(spdiags(ones(Length+1,1),0,Length,Length),H) + ...
            kron(spdiags(exp(1i*(0:Length)*varargin{2})',0,Length,Length), temp);
        H = temp + temp' + kron(Hopping_x,V) + kron(Hopping_x',V');

    elseif varargin{1} == "By"

        temp = kron(diag(exp(1i*(1:Width+1)*varargin{2})), H_intra_00);
        H = temp + temp' + kron(Hopping_y, H_intra_0y)...
            + kron(Hopping_y', H_intra_y0);
        % 由于armchair的石墨烯中，最后一层原子不存在，
        % 构建出来的哈密顿量需要删去单胞构造多出来的一层原子
        H = sparse(H(1:4*Width+2, 1:4*Width+2));
        H = kron(speye(Length), H) + kron(Hopping_x,V) + kron(Hopping_x',V');
    end
else
    % 构成单层的哈密顿量
    H = kron(eye(Width+1), H_intra_00) + kron(Hopping_y, H_intra_0y) + kron(Hopping_y', H_intra_y0);
    % 消除最后一层原子
    H = sparse(H(1:4*Width+2, 1:4*Width+2));
    % 构成层间整体哈密顿量
    H = kron(speye(Length), H) + kron(Hopping_x, V) + kron(Hopping_x', V');
end


end