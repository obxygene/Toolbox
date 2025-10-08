
function [H, H_col, T_col] = Lieb_2D_Hamiltonian(epsilon,J_plus,J_minus_x,J_minus_y,N_row,N_col, varargin)
%% LIEB_2D_HAMILTONIAN Real-space tight-binding Hamiltonian of a 2D Lieb lattice.
% H = Lieb_2D_Hamiltonian(epsilon,J_plus,J_minus_x,J_minus_y,N_row,N_col[,bc])
% epsilon   : [eA eB eC] onsite energies (1x3)
% J_plus    : intra-cell A<->B and A<->C hopping
% J_minus_x : inter-cell hopping along +x (B -> next A)
% J_minus_y : inter-cell hopping along +y (C -> next A)
% N_row     : number of unit cells in y
% N_col     : number of unit cells in x
% bc        : optional boundary condition
%             'OBC' (default), 'PBC', or struct with fields .x / .y
%             field example: bc = struct('x', 'PBC', 'y', 'OBC');
% OUTPUTS:
%   H     : full sparse Hamiltonian (3*N_row*N_col)
%   H_col : single-column block (size 3*N_row) including vertical couplings (and y-PBC if requested)
%   T_col : inter-column hopping from column j to j+1 (size 3*N_row), only B(r) -> A(r) = J_minus_x
%           (Hermitian counterpart is T_col')
%% Lieb Geometry
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...
% C   C   C   C   C ...
% A B A B A B A B A ...


% ---- checks ----
if numel(epsilon) ~= 3
    error('epsilon must be a length‑3 vector for (A,B,C) onsite energies.');
end
if ~(isscalar(J_plus) && isscalar(J_minus_x) && isscalar(J_minus_y))
    error('Hopping parameters must be scalar values.');
end
if ~(isscalar(N_row) && isscalar(N_col) && N_row>=1 && N_col>=1 && N_row==floor(N_row) && N_col==floor(N_col))
    error('N_row and N_col must be positive integers.');
end

% ---- boundary conditions ----
bc_x = "OBC"; bc_y = "OBC";

if nargin > 6 && ~isempty(varargin)
    bcArg = varargin{1};
    if (ischar(bcArg) || isstring(bcArg))
        modeStr = upper(string(bcArg));
        switch modeStr
            case "PBC"
                bc_x = "PBC"; bc_y = "PBC";
            case "OBC"
                % keep defaults
            otherwise
                error("Boundary condition string must be 'PBC' or 'OBC'.");
        end
    elseif isstruct(bcArg)
        if isfield(bcArg, 'x'); bc_x = upper(string(bcArg.x)); end
        if isfield(bcArg, 'y'); bc_y = upper(string(bcArg.y)); end
        if ~(bc_x=="OBC"||bc_x=="PBC") || ~(bc_y=="OBC"||bc_y=="PBC")
            error('Struct fields .x and .y must be ''OBC'' or ''PBC''.');
        end
    else
        error('Unsupported boundary condition specification.');
    end
end

% ---- build single-column block (H_col) and inter-column hop (T_col) if requested ----
if nargout > 1
    dcol = 3 * N_row;
    % site ordering per column: for row r (1..N_row): A= (r-1)*3+1; B=+2; C=+3
    Acol = (0:N_row-1)*3 + 1; Bcol = Acol + 1; Ccol = Acol + 2;
    % estimate nnz: onsite 3N + intra 4N + vertical 2(N_row-1) (+2N_row if PBC)
    est_col = 7*N_row + 2*(N_row-1) + (bc_y=="PBC")*2*N_row;
    Ic = zeros(est_col,1); Jc = zeros(est_col,1); Vc = zeros(est_col,1); pc = 0;
    % onsite
    Ic(pc+1:pc+N_row)=Acol; Jc(pc+1:pc+N_row)=Acol; Vc(pc+1:pc+N_row)=epsilon(1); pc=pc+N_row;
    Ic(pc+1:pc+N_row)=Bcol; Jc(pc+1:pc+N_row)=Bcol; Vc(pc+1:pc+N_row)=epsilon(2); pc=pc+N_row;
    Ic(pc+1:pc+N_row)=Ccol; Jc(pc+1:pc+N_row)=Ccol; Vc(pc+1:pc+N_row)=epsilon(3); pc=pc+N_row;
    % intra A-B & A-C
    Ic(pc+1:pc+N_row)=Acol; Jc(pc+1:pc+N_row)=Bcol; Vc(pc+1:pc+N_row)=J_plus; pc=pc+N_row;
    Ic(pc+1:pc+N_row)=Bcol; Jc(pc+1:pc+N_row)=Acol; Vc(pc+1:pc+N_row)=J_plus; pc=pc+N_row;
    Ic(pc+1:pc+N_row)=Acol; Jc(pc+1:pc+N_row)=Ccol; Vc(pc+1:pc+N_row)=J_plus; pc=pc+N_row;
    Ic(pc+1:pc+N_row)=Ccol; Jc(pc+1:pc+N_row)=Acol; Vc(pc+1:pc+N_row)=J_plus; pc=pc+N_row;
    % vertical C(r)->A(r+1)
    if N_row>1
        fromC = Ccol(1:end-1); toA = Acol(2:end); m = numel(fromC);
        Ic(pc+1:pc+m)=fromC; Jc(pc+1:pc+m)=toA; Vc(pc+1:pc+m)=J_minus_y; pc=pc+m;
        Ic(pc+1:pc+m)=toA;   Jc(pc+1:pc+m)=fromC; Vc(pc+1:pc+m)=J_minus_y; pc=pc+m;
    end
    % y-PBC wrap C(last)->A(first)
    if bc_y=="PBC" && N_row>1
        Ic(pc+1)=Ccol(end); Jc(pc+1)=Acol(1); Vc(pc+1)=J_minus_y; pc=pc+1;
        Ic(pc+1)=Acol(1);   Jc(pc+1)=Ccol(end); Vc(pc+1)=J_minus_y; pc=pc+1;
    end
    Ic=Ic(1:pc); Jc=Jc(1:pc); Vc=Vc(1:pc);
    H_col = sparse(Ic,Jc,Vc,dcol,dcol);
else
    H_col = [];
end

if nargout > 2
    if isempty(H_col) % ensure ordering exists
        dcol = 3 * N_row;
    end
    % T_col: only B(r)->A(r) entries
    Acol = (0:N_row-1)*3 + 1; Bcol = Acol + 1;
    T_col = sparse(Bcol, Acol, J_minus_x*ones(N_row,1), dcol, dcol);
else
    T_col = [];
end

% ---- vectorized full lattice assembly ----
numCells = N_row * N_col; dim = 3 * numCells;          % per cell: (A,B,C)
cells = 0:(numCells-1); A = cells*3 + 1; B = A + 1; C = A + 2;

% conservative nnz estimate
est = 12 * numCells; I = zeros(est,1); J = zeros(est,1); V = zeros(est,1); p=0;

% onsite
I(p+1:p+numCells)=A; J(p+1:p+numCells)=A; V(p+1:p+numCells)=epsilon(1); p=p+numCells;
I(p+1:p+numCells)=B; J(p+1:p+numCells)=B; V(p+1:p+numCells)=epsilon(2); p=p+numCells;
I(p+1:p+numCells)=C; J(p+1:p+numCells)=C; V(p+1:p+numCells)=epsilon(3); p=p+numCells;

% intra A-B & A-C (Hermitian)
I(p+1:p+numCells)=A; J(p+1:p+numCells)=B; V(p+1:p+numCells)=J_plus; p=p+numCells;
I(p+1:p+numCells)=B; J(p+1:p+numCells)=A; V(p+1:p+numCells)=J_plus; p=p+numCells;
I(p+1:p+numCells)=A; J(p+1:p+numCells)=C; V(p+1:p+numCells)=J_plus; p=p+numCells;
I(p+1:p+numCells)=C; J(p+1:p+numCells)=A; V(p+1:p+numCells)=J_plus; p=p+numCells;

% vertical C(i,j)->A(i+1,j)
if N_row>1
    rows=(1:(N_row-1))'; cols=1:N_col; [RR,CC]=ndgrid(rows,cols);
    from=(CC-1)*N_row+RR; to=from+1; idxC=(from-1)*3+3; idxA=(to-1)*3+1; m=numel(idxC);
    I(p+1:p+m)=idxC; J(p+1:p+m)=idxA; V(p+1:p+m)=J_minus_y; p=p+m;
    I(p+1:p+m)=idxA; J(p+1:p+m)=idxC; V(p+1:p+m)=J_minus_y; p=p+m;
end
% vertical PBC wrap last row C -> first row A
if bc_y=="PBC" && N_row>1
    last=(1:N_col)*N_row; first=last-(N_row-1); idxC=(last-1)*3+3; idxA=(first-1)*3+1; m=numel(idxC);
    I(p+1:p+m)=idxC; J(p+1:p+m)=idxA; V(p+1:p+m)=J_minus_y; p=p+m;
    I(p+1:p+m)=idxA; J(p+1:p+m)=idxC; V(p+1:p+m)=J_minus_y; p=p+m;
end

% horizontal B(i,j)->A(i,j+1)
if N_col>1
    cols=(1:(N_col-1))'; rows=1:N_row; [CC,RR]=ndgrid(cols,rows);
    from=(CC-1)*N_row+RR; to=from+N_row; idxB=(from-1)*3+2; idxA=(to-1)*3+1; m=numel(idxB);
    I(p+1:p+m)=idxB; J(p+1:p+m)=idxA; V(p+1:p+m)=J_minus_x; p=p+m;
    I(p+1:p+m)=idxA; J(p+1:p+m)=idxB; V(p+1:p+m)=J_minus_x; p=p+m;
end
% horizontal PBC wrap last col B -> first col A
if bc_x=="PBC" && N_col>1
    lastBase=(N_col-1)*N_row; lastCells=lastBase+(1:N_row); firstCells=1:N_row;
    idxB=(lastCells-1)*3+2; idxA=(firstCells-1)*3+1; m=numel(idxB);
    I(p+1:p+m)=idxB; J(p+1:p+m)=idxA; V(p+1:p+m)=J_minus_x; p=p+m;
    I(p+1:p+m)=idxA; J(p+1:p+m)=idxB; V(p+1:p+m)=J_minus_x; p=p+m;
end

% build sparse
I=I(1:p); J=J(1:p); V=V(1:p);
H = sparse(I,J,V,dim,dim);
end