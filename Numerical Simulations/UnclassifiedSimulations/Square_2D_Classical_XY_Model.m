%[text] 作者：大黄猫
%[text] 链接：[\[小练习\]二维XY模型的蒙特卡洛计算 - 知乎](https://zhuanlan.zhihu.com/p/20919617885)

Lx=20;
Ly=20;
J=1;
h=0;
% Tc=J/(2/pi)=1.5708*J
N_eq=1250*5; % number of times to equilibrium state
N_me=5000*5; % number of times to measure the state

N=39;
T=zeros(N,1);
E_density=zeros(N,1);
Cv=zeros(N,1);
susceptibility=zeros(N,1);
structure_factor=zeros(N,1);

tic
for j=1:N
    T(j)=0.1+1.9*(j-1)/(N-1);
    [E,E2,M,chi,str_f,theta_f,E_sweep]=XY2D(Lx,Ly,J,h,T(j),N_eq,N_me);
    E_density(j)=E/(Lx*Ly);              % 能量密度
    Cv(j)=((E2-E^2)/T(j)^2)/(Lx*Ly);     % 平均每格点的比热
    susceptibility(j)=chi/(Lx*Ly)^2;     % spin susceptibility
    structure_factor(j)=str_f/(Lx*Ly)^2; % spin susceptibility
end
toc

figure(1)
plot(T,E_density,'-o');
figure(2)
plot(T,Cv,'--s');
%figure(3)
%plot(T,susceptibility,'--s');
figure(4)
plot(T,structure_factor,'--s');

%TT=zeros(N-1,1);
%divS=zeros(N-1,1);
%for j=1:N-1
%    TT(j)=T(j);
%    divS(j)=(structure_factor(j+1)-structure_factor(j))/(1.9/(N-1));
%end
%figure(5)
%plot(TT,-divS,'--s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E,E2,M,chi,structure_factor_mean,theta_f,E_sweep] = XY2D(Lx,Ly,J,h,T,N_eq,N_me)
%% starting 初始化自旋构型
theta_0=starting_simulation(Lx,Ly,J,h,T); % 初始化自旋构型
% E0=energy(Lx,Ly,J,h,T,theta_0);           % 计算初始构型的能量
N_sites=Lx*Ly;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let system to equilibrium state
E_sweep=zeros(N_eq,1);     % 每一次sweep后的能量
for loop=1:N_eq            % 使得体系达到平衡需要的sweep的步数,每一个sweep需要遍历每一个格点
    for r=1:N_sites      % 遍历每一个格点
        theta=theta_0;
        [theta_out,~]=flip_site(Lx,Ly,J,h,T,theta,r);
        theta_0=theta_out;
    end
    E_sweep(loop,1)=energy(Lx,Ly,J,h,T,theta);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% measure the system, sweep one time, measure once
E_final=zeros(N_me,1);
magnet_final=zeros(N_me,1);
structure_factor_final=zeros(N_me,1);
for loop=1:N_me      % 达到平衡后,sweep的次数

    for r=1:N_sites
        theta=theta_0;
        [theta_out,magnet_out]=flip_site(Lx,Ly,J,h,T,theta,r);
        theta_0=theta_out;
        magnet_0=magnet_out;
    end
    if mod(loop,5)==1         %每间隔5个sweep计算一次能量和磁化
        E_final(loop,1)=energy(Lx,Ly,J,h,T,theta);
        magnet_final(loop,1)=magnet_0;
        structure_factor_final(loop,1)=structure_factor(Lx,Ly,J,h,T,theta);
    end
end
theta_f=theta; % manifold of XY theta

E=sum(E_final)/(N_me/5); % the mean value of energy
E2=sum(E_final.^2)/(N_me/5); % the mean value of energy^2
M=sum(magnet_final)/(N_me/5); % the mean value of magnetization
structure_factor_mean=sum(structure_factor_final)/(N_me/5);
chi=(structure_factor_mean-M^2)/T; % the mean value of susceptibility
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta] = starting_simulation(Lx,Ly,J,h,T)
%% generate 2D XY model theta(Lx*Ly,1)
N_sites=Lx*Ly;
%% generate initial manifold
theta=2*pi*rand(N_sites,1); % using random number to obtain random manifold for disordered state
% theta=2*pi*ones(N_sites,1); % theta all 2*pi for ordered state

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,magnet]=flip_site(Lx,Ly,J,h,T,theta,r)
%% Try to flip the r-th site in Ising model   简化的公式计算能量差值Delta_E,比flip_site2那个快很多!
%% output energy E and Ising manifold sigma

% calulate energy difference between fliping an Ising spin
theta_old=theta;

detaq=1; % tuning this factor to obtain good result
theta(r)=theta(r)+2*pi*detaq*rand();

theta(r)=mod(theta(r),2*pi); % 把角度变量压缩在(0,2*pi)内

[ry,rx]=site_label(r,Lx);          % r=(rx,ry)
if rx>1 && rx<Lx && ry>1 && ry<Ly
    r_right=rx+1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif ry==1 && rx>1 && rx<Lx
    r_right=rx+1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(Ly+0-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif ry==Ly && rx>1 && rx<Lx
    r_right=rx+1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(1+0-1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==1 && ry>1 && ry<Ly
    r_right=rx+1+(ry-1)*Lx;
    r_left=Lx+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==Lx && ry>1 && ry<Ly
    r_right=1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==1 && ry==1
    r_right=rx+1+(ry-1)*Lx;
    r_left=Lx+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(Ly-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==Lx && ry==1
    r_right=1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(ry-1+1)*Lx;
    r_up=rx+(Ly-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==1 && ry==Ly
    r_right=rx+1+(ry-1)*Lx;
    r_left=Lx+(ry-1)*Lx;
    r_donw=rx+(1+0-1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
elseif rx==Lx && ry==Ly
    r_right=1+(ry-1)*Lx;
    r_left=rx-1+(ry-1)*Lx;
    r_donw=rx+(1+0-1)*Lx;
    r_up=rx+(ry-1-1)*Lx;
    Delta_E=-J*(cos(theta(r)-theta_old(r_right))+cos(theta(r)-theta_old(r_left))+cos(theta(r)-theta_old(r_donw))+cos(theta(r)-theta_old(r_up)))-h*cos(theta(r))+...
        +J*(cos(theta_old(r)-theta_old(r_right))+cos(theta_old(r)-theta_old(r_left))+cos(theta_old(r)-theta_old(r_donw))+cos(theta_old(r)-theta_old(r_up)))+h*cos(theta_old(r));
end
% change or not
rate=rand(1,1);
if rate>=exp(-Delta_E/T)
    theta(r)=theta_old(r);
end
% measure magnetization
magnet=magnetization(Lx,Ly,J,h,T,theta);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=site_label(r,Lx)

% Generate the site label (a,b) for the r-th sites.
% Input:
%   Lx: The number of lattice sites in the x direction.

if rem(r,Lx)==0

    b=Lx;
    a=(r-b)/Lx+1;

else
    b=rem(r,Lx);
    a=(r-b)/Lx+1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mx=magnetization(Lx,Ly,J,h,T,theta)
%% calculate the magnetization of 2D XY model cos(theta)
Mx=sum(cos(theta));
%My=sum(sin(theta));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = energy(Lx,Ly,J,h,T,theta)
%% calculate the energy of 2D XY model theta
energy=0;

% energy of Lx-1 site

for ry=1:Ly
    for rx=1:Lx
        r=rx+(ry-1)*Lx;
        if rx<Lx && ry<Ly
            rr=rx+1+(ry-1)*Lx;
            rd=rx+(ry-1+1)*Lx;
            energy=energy-J*cos(theta(r,1)-theta(rr,1))-J*cos(theta(r,1)-theta(rd,1))-h*cos(theta(r,1));
        elseif rx<Lx && ry==Ly
            rr=rx+1+(ry-1)*Lx;
            rd=rx+(ry-1+1-Ly)*Lx;
            energy=energy-J*cos(theta(r,1)-theta(rr,1))-J*cos(theta(r,1)-theta(rd,1))-h*cos(theta(r,1));
        elseif rx==Lx && ry<Ly
            rr=rx+1-Lx+(ry-1)*Lx;
            rd=rx+(ry-1+1)*Lx;
            energy=energy-J*cos(theta(r,1)-theta(rr,1))-J*cos(theta(r,1)-theta(rd,1))-h*cos(theta(r,1));
        else
            rr=rx+1-Lx+(ry-1)*Lx;
            rd=rx+(ry-1+1-Ly)*Lx;
            energy=energy-J*cos(theta(r,1)-theta(rr,1))-J*cos(theta(r,1)-theta(rd,1))-h*cos(theta(r,1));
        end
    end
end

y=energy;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structure_factor=structure_factor(Lx,Ly,J,h,T,theta)
% structure_factor=sum_{j,k}S_j*S_k=sum_{j,k}cos(theta_j-theta_k)
structure_factor=0;
N_sites=Lx*Ly;
for j=1:N_sites
    for k=1:N_sites
        structure_factor=structure_factor+cos(theta(j)-theta(k));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
