
function Mat = get_Nearest_Square_Hopping_Sparse(col_tot, row_tot, t_x, t_y,varargin)


%%   生成二维方格子上最近邻跃迁中非零的稀疏矩阵
%    Generate the Nearest hopping nonzero sparse matrix
%   输入:
%   col_tot:方格子的列数
%   row_tot:方格子的行数
%   tx: 自左向右的跃迁积分/强度参数
%     1→tx→2 
%   ty: 自下向上的跃迁积分/强度参数
%     1
%     ⬆ ty
%    L+1
%   输出对应最近邻跃迁中不为0的稀疏矩阵
%   需要使用时可使用full函数获得完整矩阵
%   varargin{1} 不为空时给出周期性边界条件
%   varargin{2} 以朗道规范(0, Bx, 0)给出磁场
%   varargin{3} 当该参数不为零时，给出(-By, 0, 0)的磁场
% % If only 6 arguments, we choose Landau's gauge A = (0, Bx, 0)
        % If the number of arguments > 7, we choose A = (-By, 0, 0)
% 
    dim = col_tot * row_tot;
    row_list = [];
    col_list = [];
    val_list = [];
    %% Generate the hopping matrix in nearest neighbor
    % Generate index of upper part used to generate the sparse matrix
%         1   2   ...   m+1 ... col_tot
%    1    0   1    0 ... 0   
%    n    0   0   ...    1 .... 0  
%   n+1   0
%   ...   0
% row_tot 0

%% If number of argument smaller than 6, this model has no magnetic field.
    if nargin < 6
        for m = 1:row_tot
            for n = 1:col_tot-1
    
                row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                col_list = horzcat(col_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                val_list = horzcat(val_list, t_x);
    
                row_list = horzcat(row_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                val_list = horzcat(val_list, t_x');
            end
        end
    % Generate index of upper part used to generate the sparse matrix
    %         1   2   ...    m ... col_tot
    %    1    0   0    0 ... 0   
    %    2    1   0    0 ... 0   
    %    n    0   0   ...    0 ... 0  
    %   n+1   0   0   ...    1
    %   ...   0
    % row_tot 0
        for m = 1:row_tot-1
            for n = 1:col_tot
                row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                col_list = horzcat(col_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                val_list = horzcat(val_list, t_y');
    
                row_list = horzcat(row_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                val_list = horzcat(val_list, t_y);
            end
        end
    %% If there is magnetic field, add corresponding phase
    else
        % If only 6 arguments, we choose Landau's gauge A = (0, Bx, 0)
        % If the number of arguments > 7, we choose A = (-By, 0, 0)
        if nargin < 7
            for m = 1:row_tot
                for n = 1:col_tot-1
        
                    row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_x);
        
                    row_list = horzcat(row_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_x');
                end
            end

            for m = 1:row_tot-1
                for n = 1:col_tot
                    row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                    val_list = horzcat(val_list, (t_y * exp(1i * varargin{2}*(n)))');
        
                    row_list = horzcat(row_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_y * exp(1i * varargin{2}*(n)));
                end
            end
            % (By,0,0)
        else
             for m = 1:row_tot
                for n = 1:col_tot-1
        
                    row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_x * exp(1i * varargin{2}*(m)));
        
                    row_list = horzcat(row_list, siteindex2Mat(n+1, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                    val_list = horzcat(val_list, (t_x * exp(1i * varargin{2}*(m)))');
                end
            end

            for m = 1:row_tot-1
                for n = 1:col_tot
                    row_list = horzcat(row_list, siteindex2Mat(n, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                    val_list = horzcat(val_list, t_y');
        
                    row_list = horzcat(row_list, siteindex2Mat(n, m+1, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_y);
                end
            end


        end
    end




        
    %% If hopping imposed periodic boundary condition, the boundary block should be 1
        % 边界条件
        % Periodic Boundary Block
        if nargin > 4
            if ~isempty(varargin{1}) 
                for m = 1:row_tot
                    row_list = horzcat(row_list, siteindex2Mat(1, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(col_tot, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_x);
        
                    row_list = horzcat(row_list, siteindex2Mat(col_tot, m, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(1, m, col_tot, row_tot));
                    val_list = horzcat(val_list, t_x');
                end
                for n = 1:col_tot
                    row_list = horzcat(row_list, siteindex2Mat(n, 1, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, row_tot, col_tot, row_tot));
                    val_list = horzcat(val_list, t_y');
        
                    row_list = horzcat(row_list, siteindex2Mat(n, row_tot, col_tot, row_tot));
                    col_list = horzcat(col_list, siteindex2Mat(n, 1, col_tot, row_tot));
                    val_list = horzcat(val_list, t_y);
                end
            end
        end

    Mat = sparse(row_list, col_list, val_list, dim, dim);
end

%% siteindex2Mat 
function label = siteindex2Mat(col_site,row_site,col_tot,row_tot)
    label = col_site + col_tot*(row_site-1);
end