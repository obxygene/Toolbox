function [H_sparse, a, b, varargout] = KPM_scaleHamiltonian(H_sparse, varargin)
%% Kernel Polynomial Method - Scale
% Scale the Hamiltonian to make the spectrum fall in [-1, 1]
% 
% input:
% H_sparse[Matrix]: input sparse Hamiltonian
% varargin{1} epsilon[positive number]: Scale the spectrum to make sure it fall in [-1,1]
% varargin{2} ScaleFactor[positive number]: Scale the spectrum once more to avoid the
%                          boundary inaccuracy
% 
% default epsilon = 0.05 
% default ScaleFactor = 1
% 
% 
% 
%% 
    if nargin > 1
        if isempty(varargin{1})
            epsilon = 0.05;
        else
            epsilon = varargin{1};
        end
    else
        epsilon = 0.05;
    end
    if nargin > 2
        if isempty(varargin{2})
            ScaleFactor = 1;
        else
            ScaleFactor = varargin{2};
        end
    else
        ScaleFactor = 1;
    end
    if nargin > 3
        tol = epsilon / 2;
        Bounds = varargin{3};
        E_min = Bounds(1);
        E_max = Bounds(2);
    else
        tol = epsilon / 2;
        v0 = exp(2*pi*1i * rand(size(H_sparse,1), 1));
        E_max = ScaleFactor*eigs(H_sparse, 1, 'largestreal','FailureTreatment','keep','Tolerance',tol,'StartVector',v0);
        E_min = ScaleFactor*eigs(H_sparse, 1, 'smallestreal','FailureTreatment','keep','Tolerance',tol,'StartVector',v0);
    end

    a = abs(E_max - E_min)/(2 - epsilon);
    b = (E_max + E_min)/2;
    if E_max - E_min <= abs(E_max + E_min) * tol / 2
        error("The Hamiltonian has a single eigenvalue, it is not possible to obtain a spectral density.")
    end
    H_sparse = (H_sparse - b*speye(size(H_sparse)))/a;
    varargout{1} = E_min;
    varargout{2} = E_max;
    varargout{3} = epsilon;
end