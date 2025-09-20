%[text] 作者：大黄猫
%[text] 链接：[\[小练习\]Hubbard模型中的螺旋磁有序(spiral magnetic order)\[内有matlab代码\] - 知乎](https://zhuanlan.zhihu.com/p/8585399339)
%[text] 来源：知乎
%[text] 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。
%[text] [PhysRevB.81.094407](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.81.094407)
%[text] [Conditions for the Spin-Spiral State in Itinerant Magnets | Scientific.Net](https://www.scientific.net/SSP.152-153.559)

% main program 计算二维Hubbard模型在spiral态下的性质,可以参考文献[PHYSICAL REVIEW B 81,094407(2010);Solid State Phenomena Vols 152-153 (2009) pp 559-562]
clear
global t t1 t2 U T nc

t=1;
t1=0*t;
t2=0*t;
U=4.0*t;
nc=0.95;
T=0.005*t; % 温度太低会导致收敛性问题
%% 由于方程的解很多,因此具有初值敏感性,需要仔细调节初值！
%x0=[0.2;0.8;pi-0.0;pi-0.5];% for nc near 1 (pi,q)-state
x0=[0.2;0.8;pi-0.5;pi-0.5]; % for nc near 1 (q,q)-state
%x0=[0.1;-1;1;0];  % for nc near 0.3 (q,0)-state
%x0=[0.1;0.5;pi;0]; % for nc near 0.5 (pi,0)-state
options=optimset('Display','iter','MaxFunEvals',400); % Option to display output
%% 求解四个平均场方程
[x,~]=fsolve(@F1,x0,options);    % Call solver
x                                % x(1): m^2 x(2):mu x(3):Qx x(4):Qy
y=energy(x)                      % 计算能量
%% 计算磁性为零时的能量
x0=U/2;
[x2,~]=fsolve(@F2,x0,options);   % Call solver
y2=energy2(x2)                   % 计算能量
%% 通过求解平均场方程并计算不同类型磁有序态的能量就可以得到相图,例如Solid State Phenomena Vols 152-153 (2009) pp 559-562l里面的图1
%% 或者PHYSICAL REVIEW B 81,094407(2010)的图1到图4.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=DF(x,T)
% 费米部分函数
tp1=exp(x/T)+1;
y=1./tp1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=energy(A)
%% 平均场近似下体系的总能量密度 A(1):m^2 A(2):mu
global t t1 t2 U T nc
S=1/(2*pi)^2;
dd=1e-7;
    function [z1]=nf1(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A(2)+U/2*nc;
        ek2=-2*t*(cos(k1+A(3))+cos(k2+A(4)))-4*t1*cos(k1+A(3)).*cos(k2+A(4))-...
            4*t2*(cos(k1+A(3)).^2+cos(k2+A(4)).^2-1)-A(2)+U/2*nc;
        ak=sqrt((ek1-ek2).^2+U^2*A(1)+dd);
        eu=(ek1+ek2+ak)/2;
        el=(ek1+ek2-ak)/2;
        
        tp1=el.*DF(el,T)+eu.*DF(eu,T);
        
        z1=S*tp1; % 费米子部分的能量密度
    end

y=quad2d(@nf1,-pi,pi,-pi,pi)+U/4*(A(1)-nc^2); %加上常数部分的能量密度

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=energy2(A)
%% 磁性序参量为零的总能量密度 A:mu
global t t1 t2 U T nc
S=1/(2*pi)^2;
    function [z1]=nf1(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A+U/2*nc;
        
        tp1=2*ek1.*DF(ek1,T);
        
        z1=S*tp1;
    end

y=quad2d(@nf1,-pi,pi,-pi,pi);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=F1(A)
%% 求解四个平均场方程 A(1):m^2 A(2):mu A(3):Qx A(4):Qy
global t t1 t2 U T nc
S=1/(2*pi)^2;
dd=1e-7;
    function [z1]=nf1(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A(2)+U/2*nc;
        ek2=-2*t*(cos(k1+A(3))+cos(k2+A(4)))-4*t1*cos(k1+A(3)).*cos(k2+A(4))-...
            4*t2*(cos(k1+A(3)).^2+cos(k2+A(4)).^2-1)-A(2)+U/2*nc;
        ak=sqrt((ek1-ek2).^2+U^2*A(1)+dd);
        
        eu=(ek1+ek2+ak)/2;
        el=(ek1+ek2-ak)/2;
        
        tp1=(DF(el,T)-DF(eu,T))./ak;
        
        z1=S*tp1;
    end

Sp1=quad2d(@nf1,-pi,pi,-pi,pi);

y1=U*Sp1-1;

    function [z2]=nf2(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A(2)+U/2*nc;
        ek2=-2*t*(cos(k1+A(3))+cos(k2+A(4)))-4*t1*cos(k1+A(3)).*cos(k2+A(4))-...
            4*t2*(cos(k1+A(3)).^2+cos(k2+A(4)).^2-1)-A(2)+U/2*nc;
        ak=sqrt((ek1-ek2).^2+U^2*A(1)+dd);
        
        eu=(ek1+ek2+ak)/2;
        el=(ek1+ek2-ak)/2;
        
        tp1=(DF(el,T)+DF(eu,T));
        
        z2=S*tp1;
    end
Sp2=quad2d(@nf2,-pi,pi,-pi,pi);
y2=Sp2-nc;

    function [z3]=nf3(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A(2)+U/2*nc;
        ek2=-2*t*(cos(k1+A(3))+cos(k2+A(4)))-4*t1*cos(k1+A(3)).*cos(k2+A(4))-...
            4*t2*(cos(k1+A(3)).^2+cos(k2+A(4)).^2-1)-A(2)+U/2*nc;
        ak=sqrt((ek1-ek2).^2+U^2*A(1)+dd);
        
        eu=(ek1+ek2+ak)/2;
        el=(ek1+ek2-ak)/2;
        
        pek2=2*t*sin(k1+A(3))+4*t1*sin(k1+A(3)).*cos(k2+A(4))+...
            8*t2*cos(k1+A(3)).*sin(k1+A(3));
        pEu=pek2.*(1-(ek1-ek2)./ak)/2;
        pEl=pek2.*(1+(ek1-ek2)./ak)/2;
        
        tp1=(DF(el,T).*pEl+DF(eu,T).*pEu);
        
        z3=S*tp1;
    end
Sp3=quad2d(@nf3,-pi,pi,-pi,pi);
y3=Sp3;

    function [z4]=nf4(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A(2)+U/2*nc;
        ek2=-2*t*(cos(k1+A(3))+cos(k2+A(4)))-4*t1*cos(k1+A(3)).*cos(k2+A(4))-...
            4*t2*(cos(k1+A(3)).^2+cos(k2+A(4)).^2-1)-A(2)+U/2*nc;
        ak=sqrt((ek1-ek2).^2+U^2*A(1)+dd);
        
        eu=(ek1+ek2+ak)/2;
        el=(ek1+ek2-ak)/2;
        
        pek2=2*t*sin(k2+A(4))+4*t1*sin(k2+A(4)).*cos(k1+A(3))+...
            8*t2*cos(k2+A(4)).*sin(k2+A(4));
        pEu=pek2.*(1-(ek1-ek2)./ak)/2;
        pEl=pek2.*(1+(ek1-ek2)./ak)/2;
        
        tp1=(DF(el,T).*pEl+DF(eu,T).*pEu);
        
        z4=S*tp1;
    end
Sp4=quad2d(@nf4,-pi,pi,-pi,pi);
y4=Sp4;

y=[y1;y2;y3;y4];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=F2(A)
% 求解化学势 A:mu
global t t1 t2 U T nc
S=1/(2*pi)^2;
    function [z1]=nf1(k1,k2)
        ek1=-2*t*(cos(k1)+cos(k2))-4*t1*cos(k1).*cos(k2)-...
            4*t2*(cos(k1).^2+cos(k2).^2-1)-A+U/2*nc;     
        tp1=2*DF(ek1,T);
        
        z1=S*tp1;
    end
Sp1=quad2d(@nf1,-pi,pi,-pi,pi);
y1=Sp1-nc;

y=y1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
