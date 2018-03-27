
% m=0.05;
% M=0.15;
% c=0.2;
m=0.5;
M=1.5;
c=2;
r=(M-m)/(c-m);

npt_quad=4;
[qpt, qwt] = QuadLG(npt_quad);
num1=5000;
dxE=(pi/2-0)/num1;
xE=0:dxE:pi/2;
E=0;K=0;
for ne=1:num1
    for ik=1:npt_quad
        xtemp=(qpt(ik)-(-1))/2*dxE+xE(ne);
        Etemp=sqrt(1-r*(sin(xtemp))^2);
%         Ktemp=(1-r*(sin(xtemp))^2)^(3/2);
        Ktemp=(1-r*(sin(xtemp))^2)^(-1/2);
        E=E+Etemp*qwt(ik)*dxE/2;
        K=K+Ktemp*qwt(ik)*dxE/2;
    end
end
% mu_u0=2*(sqrt(2*(c-m))*(c*E-(c-M)*K))^(2/3);
% period=2*sqrt(2*(c-m)/mu_u0)*E;
mu_u0 = (2/3*sqrt(2*(c-m))*((2*(m+M)-c)*E+(c-M)*K))^(2/3);
period = 2*sqrt(2*(c-m)/mu_u0)*E;


dxE=(pi/4-0)/num1;
xE=0:dxE:pi/4;
E_mid=0;
for ne=1:num1
    for ik=1:npt_quad
        xtemp=(qpt(ik)-(-1))/2*dxE+xE(ne);
        Etemp=sqrt(1-r*(sin(xtemp))^2);
        E_mid=E_mid+Etemp*qwt(ik)*dxE/2;
    end
end
x_mid=sqrt(2*(c-m)/mu_u0)*E_mid; 
u_mid=m+(M-m)*(sin(pi/4))^2;  

Nu=10000;
dx1=(x_mid-0)/(Nu);
dx2=(period/2-x_mid)/(Nu);
deltax=min(dx1,dx2);
Xexc=[-period/2:dx2:-x_mid,-x_mid+dx1:dx1:0,0+dx1:dx1:x_mid,x_mid+dx2:dx2:period/2];
uexc=zeros(1,Nu*4+1);
rexc=zeros(1,Nu*4+1);
uexc(1,3*Nu+1)=u_mid;
uexc(1,Nu+1)=u_mid;

for n=1:Nu
    uL=uexc(1,(3*Nu+1)-(n-1));
    uR=uexc(1,(3*Nu+1)+n-1);
    
    L1=abs(sqrt(2*mu_u0*(M-uL)*(uL-m)/(c-uL)));     
    uL1=uL+(-dx1/2)*L1;
    R1=abs(sqrt(2*mu_u0*(M-uR)*(uR-m)/(c-uR)));
    uR1=uR+(dx2/2)*R1;
    
    L2=abs(sqrt(2*mu_u0*(M-uL1)*(uL1-m)/(c-uL1)));
    uL2=uL+(-dx1/2)*L2;
    R2=abs(sqrt(2*mu_u0*(M-uR1)*(uR1-m)/(c-uR1)));
    uR2=uR+(dx2/2)*R2;

    L3=abs(sqrt(2*mu_u0*(M-uL2)*(uL2-m)/(c-uL2)));
    uL3=uL+(-dx1)*L3;
    R3=abs(sqrt(2*mu_u0*(M-uR2)*(uR2-m)/(c-uR2)));
    uR3=uR+(dx2)*R3;
    
    L4=abs(sqrt(2*mu_u0*(M-uL3)*(uL3-m)/(c-uL3)));
    R4=abs(sqrt(2*mu_u0*(M-uR3)*(uR3-m)/(c-uR3)));
    
    UL=uL+(-dx1/6)*(L1+2*L2+2*L3+L4);
    UR=uR+(dx2/6)*(R1+2*R2+2*R3+R4);
    
    uexc(1,(3*Nu+1)+n)=UR;
    uexc(1,(3*Nu+1) -n)=UL;
    uexc(1,(Nu+1)+n)=UL;
    uexc(1,(Nu+1) -n)=UR;
    
end
rexc(2*Nu+1) =sqrt(abs(2*mu_u0*(M-uexc(2*Nu+1))*(uexc(2*Nu+1)-m)/(c-uexc(2*Nu+1))));
for n=1:2*Nu
    rexc(2*Nu+1+n) = sqrt(abs(2*mu_u0*(M-uexc(2*Nu+1+n))*(uexc(2*Nu+1+n)-m)/(c-uexc(2*Nu+1+n))));
    rexc(2*Nu+1-n) = -1*rexc(2*Nu+1+n);
end

% 
%  plot(Xexc,rexc,'b.-')
%  hold on
%  plot(Xexc,uexc,'r*-')
clearvars -except Xexc uexc rexc period c mu_u0
%     


