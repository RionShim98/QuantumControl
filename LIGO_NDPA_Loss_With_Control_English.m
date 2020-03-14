set(0,"DefaultTextInterpreter","latex")
set(0,"DefaultAxesFontSize",20)
set(0,"DefaultAxesFontName","Times new roman")
clear
%% Changing parameters
v=0;
Qk=eye(12);
Rk=1;
lam=3*10^6;
%% LIGO parameters
hbar=1.055*10^(-34);
L=4000;M=40;Wm=1;c=3*10^8;
ITM=2*pi*200;
lamlaser=1064*10^(-9);
Wo=2*pi*c/lamlaser;
%% Loss 
Tlos1=5*10^(-3);
Tlos2=0;
Tlos3=10^(-4);
Tlos4=10^(-3);
L1=1.5;
L2=L1;
L3=300;
Ln=0.5;
los1=c*Tlos1/L1% a_1
los2=c*Tlos2/L2% a_2
los3=c*Tlos3/L3% a_3(MCC)
los4=c*Tlos4/Ln% a_4(loop)
%% NDPA parameters;
kap1=2*c/L;n1=2.01;gam=n1*lam;
n2=0;kap2=n2*kap1;gNI=sqrt(c*gam/(2*L));
g24=sqrt(c*gam/Ln);g34=sqrt(c*kap1/Ln);% Ln‚Íƒ‹?[ƒv‚Ì’·‚³
%% LIGO power
Po=8*10^5;Rso=0;theo=0;% Original LIGO
%Pd=8*10^5;Rsd=0.1;thed=1.1;% Detuned LIGO
Pn=8*10^5;Rsn=0.01;then=pi/10;% phase cancel
%% Original LIGO dynamics
Go=(Wo/L)*sqrt(2*Po*L/(hbar*Wo*c));
Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rso)*sin(2*theo)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));
IFO=(1-Rso)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));
A2=[0 Wm 0 0;0 0 sqrt(2)*Gm 0;0 0 -IFO/2 del;sqrt(2)*Gm 0 -del -IFO/2];%[X,P,q,p]
B2=[0 0 0;0 0 sqrt(1/(hbar*M*Wm));-sqrt(IFO) 0 0;0 -sqrt(IFO) 0];%[Qin Pin F_{GW}]
C2=[0 0 0 sqrt(IFO)];%[X,P,q,p]
D2=[0 1 0];%[Qin Pin F_{GW}]
%% Original LIGO Noise
N=5*10^3;n=9;
T=N-n;
f=zeros(T,1);noise=zeros(T,1);sql=zeros(T,1);
for k=1:T
    f(k)=n+k;
    s=1i*2*pi*f(k);
    % Transfer function from [Qin Pin F_{GW}] to Pout
    G2=C2/(s*eye(4)-A2)*B2+D2;
    noise(k)=sqrt(abs(G2(1,1)/G2(1,3)/(M*L*s^2))^2/2+...
        abs(G2(1,2)/G2(1,3)/(M*L*s^2))^2/2);
    sql(k)=sqrt(abs(hbar*(Wm^2+s^2)/(M*L^2*s^4)));% SQL
end
figure;
loglog(f,noise,'b',f,sql,'k:','Linewidth',2.5)
%legend('Typical detector (LIGO)','Stardard quantum limit (SQL)')
hold on
%% Detuned LIGO dynamics
% Go=(Wo/L)*sqrt(2*Pd*L/(hbar*Wo*c));
% Gm=Go*sqrt(hbar/(M*Wm));
% del=-sqrt(Rsd)*sin(2*thed)*ITM/(1+Rsd+2*sqrt(Rsd)*cos(2*thed));
% IFO=(1-Rsd)*ITM/(1+Rsd+2*sqrt(Rsd)*cos(2*thed));
% A2=[0 Wm 0 0;0 0 sqrt(2)*Gm 0;0 0 -IFO/2 del;sqrt(2)*Gm 0 -del -IFO/2];%[X,P,q,p]
% B2=[0 0 0;0 0 sqrt(1/(hbar*M*Wm));-sqrt(IFO) 0 0;0 -sqrt(IFO) 0];%[Qin Pin F_{GW}]
% C2=[0 0 0 sqrt(IFO)];%[X,P,q,p]
% D2=[0 1 0];%[Qin Pin F_{GW}]
%% Detuned LIGO Noises
% f=zeros(T,1);noised=zeros(T,1);
% for k=1:T
%     f(k)=n+k;
%     s=1i*2*pi*f(k);
%     G2=C2/(s*eye(4)-A2)*B2+D2;
%     noised(k)=sqrt(abs(G2(1,1)/G2(1,3)/(M*L*s^2))^2/2+...
%         abs(G2(1,2)/G2(1,3)/(M*L*s^2))^2/2);
% end
% loglog(f,noised,'m:','Linewidth',2.5)
% hold on
%% NDPA+LIGO with Control Dynamics
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));
Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));

A1=[0 Wm 0 0 0 0 0 0 0 0 0 0;...
    0 0 sqrt(2)*Gm 0 0 0 0 0 0 0 0 0;...
    0 0 -IFO/2 del 0 gNI 0 0 0 0 0 0;...
    sqrt(2)*Gm 0 -del -IFO/2 -gNI 0 0 0 0 0 0 0;...
    0 0 0 gNI -los1/2 0 lam 0 0 0 0 0;...
    0 0 -gNI 0 0 -los1/2 0 -lam 0 0 0 0;...
    0 0 0 0 lam 0 -los2/2 0 0 g24 0 0;...
    0 0 0 0 0 -lam 0 -los2/2 -g24 0 0 0;...
    0 0 0 0 0 0 0 g24 -los4/2 0 0 g34;...
    0 0 0 0 0 0 -g24 0 0 -los4/2 -g34 0;...
    0 0 0 0 0 0 0 0 0 g34 -(kap2+los3)/2 0;...
    0 0 0 0 0 0 0 0 -g34 0 0 -(kap2+los3)/2];
%[X,P,qd,pd,q1,p1,q2,p2,qb,pb,qc,pc]

Bw=[0 0 0 0 0 0 0 0 0 0 0 0 0;...
    sqrt(1/(hbar*M*Wm)) 0 0 0 0 0 0 0 0 0 0 0 0;...
    0 -sqrt(IFO) 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 -sqrt(IFO) 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 -sqrt(los1) 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 -sqrt(los1) 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 -sqrt(los2) 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 -sqrt(los2) 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 -sqrt(los4) 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 -sqrt(los4) 0 0;...
    0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3) 0;...
    0 0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3)];
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Q4los P4los Q3los P3los]

Bc=[0;0;0;0;0;0;0;0;0;0;-1;1];% u is input to a_3

% use pd for GW detection
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];%[X,P,qd,pd,q1,p1,q2,p2,qb,pb,qc,pc]
D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]

% use pd for Control
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];
D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]

% use pc for Control
%C1c=[0 0 0 0 0 0 0 0 0 0 0 sqrt(kap2)];
%D1c=[0 0 0 0 1 0 0 0 0 0 0 0 0];
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]
%% LQR and Kalman filter
[F,S,E]=lqr(A1,Bc,Qk,Rk);
F;
Ql=diag([0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);
Plant=ss(A1,Bw,C1c,D1c);
% [kest,K,P]=kalman(Plant,Ql,0,0);
% K;
% K is equivalent to the following K1:
[K,~,~]=lqr(A1',C1c',Bw*Ql*Bw',D1c*Ql*D1c',Bw*Ql*D1c');
K=K';
%% Entire Controlled Dynamics
%d/dt[x e]^T=
Atot=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];%[x e]
Bwtot=[Bw;K*D1c-Bw];
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]

% ym=
Ctot=[C1m zeros(1,12)];%[x e]
Dtot=D1m;
%[F Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]
%% Entire Controlled (Stabilized) System Noise  
f=zeros(T,1);noise1=zeros(T,1);noise0=zeros(T,1);
for k=1:T
    f(k)=n+k;
    s=1i*2*pi*f(k);
    % After Stabilization
    Gwtot=Ctot/(s*eye(24)-Atot)*Bwtot+Dtot;%[F_GW Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]
    % Before Stabilization
    Gw0=C1m/(s*eye(12)-A1)*Bw+D1m;%[F_GW Qdin Pdin Qbin Pbin Q1los P1los Q2los P2los Qblos Pblos Qclos Pclos]
    %% Noise AFTER Stabilization    
    nn=0;
    for i=2:13
        nn=nn+abs(Gwtot(1,i)/Gwtot(1,1)/(M*L*s^2))^2/2;
    end
    noise1(k)=sqrt(nn);
    %% Noise BEFORE Stabilization
%     mm=0;
%     for i=2:5
%         mm=mm+abs(Gw0(1,i)/Gw0(1,1)/(M*L*s^2))^2/2;
%     end
%     noise0(k)=sqrt(mm);
end

loglog(f,noise1,'m','Linewidth',2.5)
legend('LIGO [Eq.~(41)]','SQL','With unstable filter',...
    'interpreter','latex','location','best','box','off')
% legend('LIGO ($\gamma_\textrm{IFO}=83.15$ Hz, $\Delta_d=-318$ Hz)',...
%     'SQL','With unstable filter',...
%     'interpreter','latex')
ylim([10^(-25) 10^(-22)])
xlim([n+1 N]);
xlabel('$\Omega/2\pi$   (Hz)');
ylabel('$\sqrt{S(\Omega)}$   $(\textrm{Hz}^{-1/2})$');
%title({'Both blue and pink lines are';'$\gamma_\textrm{IFO}=1062$ Hz, $\Delta_d=-63.0$ Hz'})

Aeig=eig(A1);
BFeig=eig(A1-Bc*F);
KCeig=eig(A1-K*C1c);
Atoteig=eig(Atot)

det(ctrb(A1,Bc));
det(obsv(A1,C1c));