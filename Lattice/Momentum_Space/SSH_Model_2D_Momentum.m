function H = SSH_Model_2D_Momentum(kx,ky,epsilon,w,v)
%SSH_MODEL_2D_HAMILTONIAN 此处显示有关此函数的摘要
%   此处显示详细说明
    H_intra = [epsilon(1), w, w, 0;
               w, epsilon(2), 0, w;
               w, 0, epsilon(3), w;
               0, w, w, epsilon(4)];
    H_hori = [0,0,0,0;v,0,0,0;0,0,0,0;0,0,v,0];
    H_vect = [0,0,0,0;0,0,0,0;v,0,0,0;0,v,0,0];
    H = H_intra + exp(1i * kx) * H_hori + exp(-1i * kx) * H_hori'...
            +exp(1i * ky) * H_vect + exp(-1i * ky) * H_vect';

end

