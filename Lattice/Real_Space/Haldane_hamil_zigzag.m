function [H00,H10,H] = Haldane_hamil_zigzag(t,M,t2,phi, Ny, Nx)

    [H_00, Hy_0, Hx_01, Hx_02] = Haldane_Unitcell_zigzag(t,M,t2,phi);

    H00 = kron(speye(Ny),H_00)+...
        kron(spdiags(ones(Ny-1,1),-1,Ny, Ny)',Hy_0')+...
        kron(spdiags(ones(Ny-1,1),-1,Ny, Ny),Hy_0);
    
    H10=kron(speye(Ny),Hx_01)+...
        kron(spdiags(ones(Ny-1,1),-1,Ny, Ny)',Hx_02);
    
    H=kron(speye(Nx),H00)+...
        kron(spdiags(ones(Nx-1,1),-1,Nx, Nx)',H10')+...
        kron(spdiags(ones(Nx-1,1),-1,Nx, Nx),H10);

% 单胞跃迁函数


    function [H_00, Hy_0, Hx_01, Hx_02] = Haldane_Unitcell_zigzag(t,M,t2,phi)
    
        H_00=[M, t;
              t, -M];
        Hx_01=[t2 * exp(-1i*phi), t;
              0, t2*exp(1i*phi)];
        Hx_02=[t2 * exp(1i*phi), 0;
              0, t2 * exp(-1i * phi)];
        Hy_0=[t2 * exp(1i*phi), t;
              0, t2 * exp(-1i * phi)];
    end

end




% 作者：邵锴
% 链接：https://zhuanlan.zhihu.com/p/6107524977
% 来源：知乎
% 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。