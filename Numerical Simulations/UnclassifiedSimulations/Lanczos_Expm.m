% The Following code comes from 
% https://github.com/juanjosegarciaripoll
% https://zhuanlan.zhihu.com/p/1919645746040399200
% 
clear;
N = 10000;
A = randn(N);A = (A + A') / 2;
v = randn(N,1);
tic;
result = expm_lanczos(A, v, 'dt', 0.1, 'order', 100, 'tol', 1e-12);
toc; %[output:4158369f]
tic;
F = expm(1i*0.1*A)*v;
toc; %[output:7c9f48c5]
all(isapprox(F, result,"tight")) %[output:49dbc22e]

function output = lanczos_expm_apply(H, v, dt, order, tol)
    % Apply the Lanczos approximation of the exponential exp(1i*dt*H)
    % onto the vector v.
    
    if nargin < 3, dt = 1.0; end
    if nargin < 4, order = []; end
    if nargin < 5, tol = 1e-14; end
    
    % Determine if H is a function handle or matrix
    if isa(H, 'function_handle')
        H_func = H;
    else
        H_func = @(x) H * x;
    end
    
    % Initialize parameters exactly as in Python
    if isempty(order)
        nmax = 10;
    else
        nmax = order;
    end
    
    v = v(:);  % Ensure column vector
    if nmax > length(v)
        nmax = length(v);
    end
    
    % Construct projected version of matrix H on Krylov subspace
    v = v;  % Keep original v
    vnrm = norm(v);
    vn = v / vnrm;
    vnm1 = zeros(size(v));
    alpha = [];
    beta = [0.0];
    start_n = 1;  % MATLAB 1-based indexing
    lasterr = inf;
    
    while true
        % Iteratively extend the Krylov basis using Lanczos recursion
        n = start_n;
        while n <= nmax
            w = H_func(vn);
            alpha_n = vn' * w;  % Complex conjugate dot product (vdot)
            alpha = [alpha, alpha_n];
            
            % Use proper indexing: alpha has been appended, beta uses original indexing
            current_alpha_idx = length(alpha);
            current_beta_idx = length(beta);
            
            w = w - alpha(current_alpha_idx) * vn - beta(current_beta_idx) * vnm1;
            vnm1 = vn;
            aux = norm(w);
            beta = [beta, aux];
            
            if aux < 1e-20
                break;
            end
            vn = w / aux;
            n = n + 1;
        end
        
        % Diagonalize the banded matrix formed by alpha and beta
        % Python: w, u = scipy.linalg.eig_banded(numpy.array([β[:-1],α]))
        n_alpha = length(alpha);
        
        if n_alpha == 1
            w_eig = alpha(1);
            u = 1;
        else
            % Create banded matrix format: [sub_diagonal; main_diagonal]
            beta_sub = beta(1:end-1);  % β[:-1] in Python
            
            % In Python eig_banded, first row is sub/super diagonal, second is main
            % We need to construct the tridiagonal matrix properly
            T = diag(alpha);  % Main diagonal
            if length(beta_sub) > 1
                off_diag = beta_sub(2:end);  % Skip the first 0.0
                if ~isempty(off_diag) && length(off_diag) >= n_alpha - 1
                    diag_elements = off_diag(1:n_alpha-1);
                    T = T + diag(diag_elements, 1) + diag(diag_elements, -1);
                end
            end
            [u, D] = eig(T);
            w_eig = diag(D);
        end
        
        % Compute exponential: fHt = u @ (numpy.exp(1j*dt*w) * u[0,:].conj())
        if size(u, 1) == 1
            fHt = u * exp(1i * dt * w_eig) * conj(u);
        else
            first_row_conj = conj(u(1, :)).';  % u[0,:].conj() in Python
            exp_diag = exp(1i * dt * w_eig);
            fHt = u * (exp_diag .* first_row_conj);
        end
        
        % Estimate error: err = abs(fHt[n]*β[n+1])
        % Note: n here refers to the last iteration index
        n_last = length(alpha);
        if n_last + 1 <= length(beta)
            err = abs(fHt(n_last) * beta(n_last + 1));
        else
            err = 0;
        end
        
        if err < tol
            break;
        end
        
        if lasterr < err || nmax == length(v) || ~isempty(order)
            warning('Lanczos failed to converge at %d iterations with error %e', ...
                    length(alpha), err);
            lasterr = err;
            break;
        end
        
        start_n = nmax + 1;
        nmax = min(floor(1.5 * nmax + 1), length(v));
    end
    
    % Recompute the Lanczos basis and apply exponential
    vnm1 = v;  % Reset to original v
    vn = v;    % Reset to original v
    output = fHt(1) * vn;
    
    for n = 2:length(fHt)
        w = H_func(vn) - alpha(n-1) * vn - beta(n-1) * vnm1;
        vnm1 = vn;
        vn = w / beta(n);
        output = output + fHt(n) * vn;
    end
end

function output = expm_lanczos(A, v, varargin)
    % Apply the Lanczos approximation of the exponential exp(1i*dt*A)
    % onto the vector or matrix v.
    
    % Parse input arguments
    p = inputParser;
    addRequired(p, 'A');
    addRequired(p, 'v');
    addParameter(p, 'dt', 1.0);
    addParameter(p, 'order', []);
    addParameter(p, 'tol', 1e-14);
    parse(p, A, v, varargin{:});
    
    dt = p.Results.dt;
    order = p.Results.order;
    tol = p.Results.tol;
    
    % Handle matrix input for v
    if size(v, 2) > 1
        output = zeros(size(v));
        for i = 1:size(v, 2)
            output(:, i) = lanczos_expm_apply(A, v(:, i), dt, order, tol);
        end
    else
        output = lanczos_expm_apply(A, v, dt, order, tol);
    end
end

% Improved test function with better convergence
function test_convergence()
    % Test with different matrix types to check convergence
    n = 30;  % Smaller size for testing
    
    fprintf('Testing Lanczos convergence...\n');
    
    % Test 1: Symmetric matrix (should converge well)
    fprintf('\nTest 1: Symmetric matrix\n');
    A1 = randn(n, n);
    A1 = (A1 + A1') / 2;
    A1 = A1 / norm(A1) * 0.5;  % Scale to improve conditioning
    v1 = randn(n, 1);
    
    try
        result1 = expm_lanczos(A1, v1, 'dt', 0.1, 'order', 25, 'tol', 1e-10);
        fprintf('Symmetric test: SUCCESS, norm = %.6f\n', norm(result1));
    catch ME
        fprintf('Symmetric test: FAILED - %s\n', ME.message);
    end
    
    % Test 2: Antisymmetric matrix (often better for exp(i*dt*A))
    fprintf('\nTest 2: Antisymmetric matrix\n');
    A2 = randn(n, n);
    A2 = (A2 - A2') / 2;  % Antisymmetric
    A2 = A2 / norm(A2) * 0.5;
    v2 = randn(n, 1);
    
    try
        result2 = expm_lanczos(A2, v2, 'dt', 0.1, 'order', 25, 'tol', 1e-10);
        fprintf('Antisymmetric test: SUCCESS, norm = %.6f\n', norm(result2));
    catch ME
        fprintf('Antisymmetric test: FAILED - %s\n', ME.message);
    end
    
    % Test 3: Compare with MATLAB's expm
    fprintf('\nTest 3: Accuracy comparison\n');
    if n <= 50
        exact = expm(1i * 0.1 * A1) * v1;
        rel_err = norm(result1 - exact) / norm(exact);
        fprintf('Relative error vs expm: %.2e\n', rel_err);
    end
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:4158369f]
%   data: {"dataType":"text","outputData":{"text":"Elapsed time is 2.939851 seconds.\n","truncated":false}}
%---
%[output:7c9f48c5]
%   data: {"dataType":"text","outputData":{"text":"Elapsed time is 364.519233 seconds.\n","truncated":false}}
%---
%[output:49dbc22e]
%   data: {"dataType":"textualVariable","outputData":{"header":"logical","name":"ans","value":"   1\n"}}
%---
