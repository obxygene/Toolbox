function [X, Y] = Graphene_ZigZag_Real_Coords(Length, Width, varargin)
% Define lattice parameters
if nargin < 3
    % a = 1.42; % Carbon-carbon bond length in angstroms
    a = 1;
else
    a = varargin{1};
end

% 单胞内4个原子的相对坐标
cell_x = [0, sqrt(3)/2, sqrt(3)/2, 0]*a;
cell_y = [0, 0.5, 1.5, 2]*a;

% 单胞间距
dx = sqrt(3)*a;
dy = 3*a;

% 预分配
X = zeros((4*Width+2) * Length,1);
Y = zeros((4*Width+2) * Length,1);
idx = 1;

for n = 0:Length-1
    for m = 0:Width-1
        for k = 1:4
            X(idx) = n*dx + cell_x(k);
            Y(idx) = m*dy + cell_y(k);
            idx = idx + 1;
        end
    end
end
% 绘图（可选）
% figure; hold on; axis equal;
% scatter(X, Y, 'ko','filled');
% for i = 1:length(X)
%     text(X(i), Y(i), num2str(i), 'FontSize', 8, 'Color', 'b');
% end
% set(gca,"YDir", "reverse")
% xlabel('x'); ylabel('y'); title('Armchair格点（4原子单胞）编号');

end