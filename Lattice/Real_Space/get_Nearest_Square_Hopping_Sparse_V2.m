function Mat = get_Nearest_Square_Hopping_Sparse_V2(col_tot, row_tot, t_x, t_y,varargin)
% 生成二维方格子最近邻跃迁的稀疏矩阵（支持磁场和周期性边界）
% 输入参数:
%   col_tot: 格子的列数
%   row_tot: 格子的行数
%   t_x: 横向（x方向）跃迁强度（可为复数）
%   t_y: 纵向（y方向）跃迁强度（可为复数）
%   varargin{1}: 非空时启用周期性边界条件（PBC）
%   varargin{2}: 字符，当有输入时，应为"Bx"或"-By"
%   varargin{3}: 数值，表示磁场相位系数
% 输出:
%   Mat: 稀疏哈密顿矩阵（dim x dim），dim = col_tot * row_tot

dim = col_tot * row_tot;
max_nz = 4*col_tot*row_tot + 2*(col_tot+row_tot); % 预估最大非零元数
row_list = zeros(max_nz,1);
col_list = zeros(max_nz,1);
val_list = zeros(max_nz,1);
nz = 0;

if nargin < 6
    % 无磁场情形
    % 横向最近邻（x方向）
    for m = 1:row_tot
        for n = 1:col_tot-1
            idx1 = n + col_tot*(m-1);    % 当前格点
            idx2 = n+1 + col_tot*(m-1);  % 右侧格点
            nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_x;
            nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_x';
        end
    end
    % 纵向
    for m = 1:row_tot-1
        for n = 1:col_tot
            idx1 = n + col_tot*(m-1);
            idx2 = n + col_tot*m;
            nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_y';
            nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_y;
        end
    end
else
    % 有磁场情形
    % 判断磁场类型
    if varargin{2} == "Bx"
        B = varargin{3};
        % Landau规范(0,Bx,0)
        for m = 1:row_tot
            for n = 1:col_tot-1
                idx1 = n + col_tot*(m-1);
                idx2 = n+1 + col_tot*(m-1);
                nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_x;
                nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_x';
            end
        end
        for m = 1:row_tot-1
            for n = 1:col_tot
                idx1 = n + col_tot*(m-1);
                idx2 = n + col_tot*m;
                phase = exp(1i * B*(n));
                nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=(t_y*phase)';
                nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_y*phase;
            end
        end
    elseif varargin{2} == "-By"
        B = varargin{3};
        % (-By,0,0)
        for m = 1:row_tot
            for n = 1:col_tot-1
                idx1 = n + col_tot*(m-1);
                idx2 = n+1 + col_tot*(m-1);
                phase = exp(1i * B*(m));
                nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_x*phase;
                nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=(t_x*phase)';
            end
        end
        for m = 1:row_tot-1
            for n = 1:col_tot
                idx1 = n + col_tot*(m-1);
                idx2 = n + col_tot*m;
                nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_y';
                nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_y;
            end
        end
    end
end

% 周期性边界
if nargin > 4 && ~isempty(varargin{1})
    for m = 1:row_tot
        idx1 = 1 + col_tot*(m-1);
        idx2 = col_tot + col_tot*(m-1);
        nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_x;
        nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_x';
    end
    for n = 1:col_tot
        idx1 = n;
        idx2 = n + col_tot*(row_tot-1);
        nz=nz+1; row_list(nz)=idx1; col_list(nz)=idx2; val_list(nz)=t_y';
        nz=nz+1; row_list(nz)=idx2; col_list(nz)=idx1; val_list(nz)=t_y;
    end
end

row_list = row_list(1:nz);
col_list = col_list(1:nz);
val_list = val_list(1:nz);

Mat = sparse(row_list, col_list, val_list, dim, dim);
end