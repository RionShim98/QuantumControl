set(0,"DefaultTextInterpreter","latex")
% set(0,"DefaultAxesFontSize",28)
% set(0,"DefaultAxesFontName","Times new roman")
%% Changing parameters
v=0;
Qk=eye(12);
Rk=1;
lam=3*10^6;
%% GW detector parameters
hbar=1.055*10^(-34);
L=4000;M=40;Wm=1;c=3*10^8;
ITM=2*pi*200;
lamlaser=1064*10^(-9);
Wo=2*pi*c/lamlaser;
%% consistent NDPA paramaters 
kap1=2*c/L;n1=2.01;gam=n1*lam;
n2=0;kap2=n2*kap1;gNI=sqrt(c*gam/(2*L));
%% consistent Loss setup
Tlos1=5*10^(-3);
Tlos2=0;
Tlos3=10^(-4);
Tlos4=10^(-3);
L1=[1.5 5 10 50];
L2=L1;
L3=300;
Ln=0.5;
leg={'LIGO [Eq.~(41)]','SQL',...
    '$L_\textrm{NDPA}=1.5$ m',...
    '$L_\textrm{NDPA}=5$ m',...
    '$L_\textrm{NDPA}=10$ m',...
    '$L_\textrm{NDPA}=50$ m'};
los3=c*Tlos3/L3;% a_3 (MCC) loss
los4=c*Tlos4/Ln;% a_4 (loop) loss
g24=sqrt(c*gam/Ln);g34=sqrt(c*kap1/Ln);
%% Power & Reflectivity & Detuning
Po=8*10^5;Rso=0;theo=0;% Original LIGO
Pn=8*10^5;Rsn=0.01;then=pi/10;% LIGO with unstable filter setup
%% Original detector dynamics & noise
Go=(Wo/L)*sqrt(2*Po*L/(hbar*Wo*c));
Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rso)*sin(2*theo)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));
IFO=(1-Rso)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));

A2=[0 Wm 0 0;0 0 sqrt(2)*Gm 0;0 0 -IFO/2 del;sqrt(2)*Gm 0 -del -IFO/2];%[X,P,q,p]
B2=[0 0 0;0 0 sqrt(1/(hbar*M*Wm));-sqrt(IFO) 0 0;0 -sqrt(IFO) 0];%[Qin Pin F_{GW}]
C2=[0 0 0 sqrt(IFO)];%[X,P,q,p]
D2=[0 1 0];%[Qin Pin F_{GW}]

N=5*10^3;n=9;
T=N-n;
f=zeros(T,1);noise=zeros(T,1);sql=zeros(T,1);
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);
G2=C2/(s*eye(4)-A2)*B2+D2;
%G(1,1) is shot noise, G(1,2) is back action noise
noise(k)=sqrt(abs(G2(1,1)/G2(1,3)/(M*L*s^2))^2/2+...
    abs(G2(1,2)/G2(1,3)/(M*L*s^2))^2/2);
sql(k)=sqrt(abs(hbar*(Wm^2+s^2)/(M*L^2*s^4)));% SQL
end
figure;
loglog(f,noise,'b','Linewidth',2.5)
hold on
loglog(f,sql,'k','Linewidth',2.5)
hold on

%% loop length change simulation
Atot=nan(24,24,length(L1));
col={'b:','g--','c:','r-.'};
for i=1:length(L1)
los1=c*Tlos1/L1(i);% a_1 loss
fprintf('los1 = %d\n',los1)
los2=c*Tlos2/L2(i);% a_1 loss
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
A1=[0 Wm 0 0 0 0 0 0 0 0 0 0;0 0 sqrt(2)*Gm 0 0 0 0 0 0 0 0 0;...
    0 0 -IFO/2 del 0 gNI 0 0 0 0 0 0;sqrt(2)*Gm 0 -del -IFO/2 -gNI 0 0 0 0 0 0 0;...
    0 0 0 gNI -los1/2 0 lam 0 0 0 0 0;0 0 -gNI 0 0 -los1/2 0 -lam 0 0 0 0;...
    0 0 0 0 lam 0 -los2/2 0 0 g24 0 0;0 0 0 0 0 -lam 0 -los2/2 -g24 0 0 0;...
    0 0 0 0 0 0 0 g24 -los4/2 0 0 g34;0 0 0 0 0 0 -g24 0 0 -los4/2 -g34 0;...
    0 0 0 0 0 0 0 0 0 g34 -(kap2+los3)/2 0;0 0 0 0 0 0 0 0 -g34 0 0 -(kap2+los3)/2];
Bw=[0 0 0 0 0 0 0 0 0 0 0 0 0;sqrt(1/(hbar*M*Wm)) 0 0 0 0 0 0 0 0 0 0 0 0;...
    0 -sqrt(IFO) 0 0 0 0 0 0 0 0 0 0 0;0 0 -sqrt(IFO) 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 -sqrt(los1) 0 0 0 0 0 0 0;0 0 0 0 0 0 -sqrt(los1) 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 -sqrt(los2) 0 0 0 0 0;0 0 0 0 0 0 0 0 -sqrt(los2) 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 -sqrt(los4) 0 0 0;0 0 0 0 0 0 0 0 0 0 -sqrt(los4) 0 0;...
    0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3) 0;0 0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3)];
Bc=[0;0;0;0;0;0;0;0;0;0;-1;1];
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
[F,~,~]=lqr(A1,Bc,Qk,Rk);
Ql=diag([v^2 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);
% Plant=ss(A1,Bw,C1c,D1c);
% [kest,K,P]=kalman(Plant,Ql,0,0);
[K,~,~]=lqr(A1',C1c',Bw*Ql*Bw',D1c*Ql*D1c',Bw*Ql*D1c');
K=K';
Atot(:,:,i)=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];
Bwtot=[Bw;K*D1c-Bw];Ctot=[C1m zeros(1,12)];Dtot=D1m;
noise1=zeros(T,1);
for k=1:T
s=1i*2*pi*f(k);
Gwtot=Ctot/(s*eye(24)-Atot(:,:,i))*Bwtot+Dtot;%% w?¨ym ‚Ì“`’BŠÖ?”[F_{GW} Qdin Pdin Qbin Pbin]
nn=0;
    for Z=1:12
        nn=nn+abs(Gwtot(1,Z+1)/Gwtot(1,1)/(M*L*s^2))^2/2;
    end
noise1(k)=sqrt(nn);
end
loglog(f,noise1,col{i},'Linewidth',2.5)
hold on
end
hold off
%% Plot status
legend(leg{:},'interpreter','latex','FontSize',20,'Box','off')
ylim([10^(-25) 10^(-22)])
xlim([n+1 N]);
ax=gca;
ax.FontSize=28;
ax.FontName='Times new roman';
xlabel('$\Omega/2\pi$   (Hz)');
ylabel('$\sqrt{S(\Omega)}$   (Hz$^{-1/2}$)');
title('(a)')
%% Staility
eigAtot=nan(24,length(L1));
for i=1:length(L1)
    eigAtot(:,i)=eig(Atot(:,:,i));
    stab=real(eigAtot(:,i))>0;
    stapoles=length(stab(stab==false));
    if stapoles==length(stab)
        fprintf('%d Stable \n',i);
    else
        fprintf('%d Unstable, %d unstable poles\n',i,length(stab)-stapoles);
    end
end
