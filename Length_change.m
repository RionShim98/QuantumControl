%%%%’²?®‰Â”\‚È‚»‚ê‚¼‚ê‚Ìƒpƒ‰ƒ??[ƒ^%%%%%%%%%%%%%%%%%
%v=40*4000*(1/8)*10^(-21);%—ÍF_{GW}‚Ì•ªŽU‚ÌŒ©?Ï‚è
v=0;
Qk=10*eye(12);Rk=0.001;%LQR‚ÌƒRƒXƒgŠÖ?”‚Ìƒpƒ‰ƒ??[ƒ^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%

%%?d—Í”gŒŸ?oŠí‚Ìƒpƒ‰ƒ??[ƒ^
hbar=1.055*10^(-34);
L=4000;M=40;Wm=1;
ITM=2*pi*200;c=3*10^8;
lamlaser=1064*10^(-9);%ƒŒ?[ƒU?[Žü”g?”
Wo=2*pi*c/lamlaser;

%%%%% ƒ?ƒX‚Ìƒpƒ‰ƒ??[ƒ^ %%%%%%
Tlos1=5*10^(-3);
Tlos2=0;
Tlos3=10^(-4);
Tlos4=10^(-3);
L1=1.5;
L2=L1;
L3=300;
Ln=0.5;
los1=c*Tlos1/L1% a_1‚Ìƒ?ƒXƒŒ?[ƒg
los2=c*Tlos2/L2;% a_2‚Ìƒ?ƒXƒŒ?[ƒg
los3=c*Tlos3/L3;% a_3(MCC)‚Ìƒ?ƒXƒŒ?[ƒg
los4=c*Tlos4/Ln;% a_4(ƒ‹?[ƒvƒLƒƒƒrƒeƒB)‚Ìƒ?ƒXƒŒ?[ƒg
%%MCCFBNDPA‚Ìƒpƒ‰ƒ??[ƒ^
% lam‚ð•Ï‚¦‚é‚Æ–Ê”’‚¢‚±‚Æ‚ª‹N‚±‚é?B—á‚¦‚Îlam=10^5
kap1=2*c/L;lam=3*10^6;n1=2.01;gam=n1*lam;
n2=0;kap2=n2*kap1;gNI=sqrt(c*gam/(2*L));
g24=sqrt(c*gam/Ln);g34=sqrt(c*kap1/Ln);% Ln‚Íƒ‹?[ƒv‚Ì’·‚³
% Original LIGO
Po=8*10^5;Rso=0;theo=0;
% FBamp‚ÆˆÀ’è‰»?§Œä‚ð‚·‚é?Û‚ÌLIGO
Pn=8*10^5;Rsn=0.01;then=pi/10;

% Original LIGO‚Ìƒpƒ‰ƒ??[ƒ^
Go=(Wo/L)*sqrt(2*Po*L/(hbar*Wo*c));
Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rso)*sin(2*theo)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));
IFO=(1-Rso)*ITM/(1+Rso+2*sqrt(Rso)*cos(2*theo));

%%% ƒIƒŠƒWƒiƒ‹‚ÌLIGO
% -M*Wm^2‚Íƒ[ƒ?‚Æ‚Ý‚È‚µ‚Ä‚¢‚é ŽÀ?Û‚Í-Wm‚Ì?€‚ð—Ž‚Æ‚µ‚Ä‚¢‚é
A2=[0 Wm 0 0;0 0 sqrt(2)*Gm 0;0 0 -IFO/2 del;sqrt(2)*Gm 0 -del -IFO/2];%[X,P,q,p]
B2=[0 0 0;0 0 sqrt(1/(hbar*M*Wm));-sqrt(IFO) 0 0;0 -sqrt(IFO) 0];%[Qin Pin F_{GW}]
C2=[0 0 0 sqrt(IFO)];%[X,P,q,p]
D2=[0 1 0];%[Qin Pin F_{GW}]

%%?Å?IŽü”g?”N‚ÆŠJŽnŽü”g?”n
N=5*10^3;n=9;
T=N-n;
%%ˆÈ?ã‚É‚æ‚Á‚Än Hz?`N Hz‚Ìƒvƒ?ƒbƒg‚ª?o—ˆ‚é
%–Ú?·Š„‚è“–‚Ä
f=zeros(T,1);noise=zeros(T,1);sql=zeros(T,1);
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);

G2=C2/(s*eye(4)-A2)*B2+D2;%%NDPA–³‚Ì“`’BŠÖ?”

%•?’Ê‚ÌLIGO‚ÌƒmƒCƒY‚Ì‘å‚«‚³
%G(1,1)‚ªshot noise, G(1,2)‚ªback action noise
noise(k)=sqrt(abs(G2(1,1)/G2(1,3)/(M*L*s^2))^2/2+...
    abs(G2(1,2)/G2(1,3)/(M*L*s^2))^2/2);

sql(k)=sqrt(abs(hbar*(Wm^2+s^2)/(M*L^2*s^4)));% SQL
end
figure;
loglog(f,noise,'b:','Displayname','Original detector','Linewidth',2.5)
hold on
loglog(f,sql,'k:','Displayname','SQL','Linewidth',2.5)
hold on

%%%%%%%%%%%%%%%%1‰ñ–Ú%%%%%%%%%%%%%%%%%%
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
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
Bc=[1;0;0;0;0;0;0;0;0;0;0;0];
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
[F,S,E]=lqr(A1,Bc,Qk,Rk);
Ql=diag([v^2 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);Plant=ss(A1,Bw,C1c,D1c);[kest,K,P]=kalman(Plant,Ql,0,0);
Atot=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];Bwtot=[Bw;K*D1c-Bw];Ctot=[C1m zeros(1,12)];Dtot=D1m;
f=zeros(T,1);noise1=zeros(T,1);
%%% LQG?§Œä‚Ì‚Ý
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);
Gwtot=Ctot/(s*eye(24)-Atot)*Bwtot+Dtot;%% w?¨ym ‚Ì“`’BŠÖ?”[F_{GW} Qdin Pdin Qbin Pbin]
nn=0;
for Z=1:12
    nn=nn+abs(Gwtot(1,Z+1)/Gwtot(1,1)/(M*L*s^2))^2/2;
end
noise1(k)=sqrt(nn);
end
loglog(f,noise1,'Displayname','L_{loop}=0.5 m','Linewidth',2.5)
hold on

%%%%%%%%%%%%%%%%%%%%2‰ñ–Ú%%%%%%%%%%%%%%%%
Ln=5;los4=c*Tlos4/Ln
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
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
    0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3) 0;0 0 0 0 -sqrt(kap2) 0 0 0 0 0 0 0 -sqrt(los3)];Bc=[1;0;0;0;0;0;0;0;0;0;0;0];
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
[F,S,E]=lqr(A1,Bc,Qk,Rk);
Ql=diag([v^2 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);Plant=ss(A1,Bw,C1c,D1c);[kest,K,P]=kalman(Plant,Ql,0,0);
Atot=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];Bwtot=[Bw;K*D1c-Bw];Ctot=[C1m zeros(1,12)];Dtot=D1m;
f=zeros(T,1);noise1=zeros(T,1);
%%% LQG?§Œä‚Ì‚Ý
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);
Gwtot=Ctot/(s*eye(24)-Atot)*Bwtot+Dtot;%% w?¨ym ‚Ì“`’BŠÖ?”[F_{GW} Qdin Pdin Qbin Pbin]
nn=0;
for Z=1:12
    nn=nn+abs(Gwtot(1,Z+1)/Gwtot(1,1)/(M*L*s^2))^2/2;
end
noise1(k)=sqrt(nn);
end
loglog(f,noise1,'Displayname','L_{loop}=5 m','Linewidth',2.5)
hold on

%%%%%%%%%%%%%%%3‰ñ–Ú%%%%%%%%%%%%%%%%%%%%
Ln=50;los4=c*Tlos4/Ln
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
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
Bc=[1;0;0;0;0;0;0;0;0;0;0;0];
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
[F,S,E]=lqr(A1,Bc,Qk,Rk);
Ql=diag([v^2 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);Plant=ss(A1,Bw,C1c,D1c);[kest,K,P]=kalman(Plant,Ql,0,0);
Atot=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];Bwtot=[Bw;K*D1c-Bw];Ctot=[C1m zeros(1,12)];Dtot=D1m;
f=zeros(T,1);noise1=zeros(T,1);
%%% LQG?§Œä‚Ì‚Ý
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);
Gwtot=Ctot/(s*eye(24)-Atot)*Bwtot+Dtot;%% w?¨ym ‚Ì“`’BŠÖ?”[F_{GW} Qdin Pdin Qbin Pbin]
nn=0;
for Z=1:12
    nn=nn+abs(Gwtot(1,Z+1)/Gwtot(1,1)/(M*L*s^2))^2/2;
end
noise1(k)=sqrt(nn);
end
loglog(f,noise1,'Displayname','L_{loop}=50 m','Linewidth',2.5)
hold on

%%%%%%%%%%%%%%%%%4‰ñ–Ú%%%%%%%%%%%%%%%%%%%%
Ln=100;los4=c*Tlos4/Ln
Go=(Wo/L)*sqrt(2*Pn*L/(hbar*Wo*c));Gm=Go*sqrt(hbar/(M*Wm));
del=-sqrt(Rsn)*sin(2*then)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));IFO=(1-Rsn)*ITM/(1+Rsn+2*sqrt(Rsn)*cos(2*then));
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
Bc=[1;0;0;0;0;0;0;0;0;0;0;0];
C1m=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1m=[0 0 1 0 0 0 0 0 0 0 0 0 0];
C1c=[0 0 0 sqrt(IFO) 0 0 0 0 0 0 0 0];D1c=[0 0 1 0 0 0 0 0 0 0 0 0 0];
[F,S,E]=lqr(A1,Bc,Qk,Rk);
Ql=diag([v^2 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]);Plant=ss(A1,Bw,C1c,D1c);[kest,K,P]=kalman(Plant,Ql,0,0);
Atot=[A1-Bc*F -Bc*F;zeros(12) A1-K*C1c];Bwtot=[Bw;K*D1c-Bw];Ctot=[C1m zeros(1,12)];Dtot=D1m;
f=zeros(T,1);noise1=zeros(T,1);
%%% LQG?§Œä‚Ì‚Ý
for k=1:T
f(k)=n+k;
s=1i*2*pi*f(k);
Gwtot=Ctot/(s*eye(24)-Atot)*Bwtot+Dtot;%% w?¨ym ‚Ì“`’BŠÖ?”[F_{GW} Qdin Pdin Qbin Pbin]
nn=0;
for Z=1:12
    nn=nn+abs(Gwtot(1,Z+1)/Gwtot(1,1)/(M*L*s^2))^2/2;
end
noise1(k)=sqrt(nn);
end
loglog(f,noise1,'Displayname','L_{loop}=100 m','Linewidth',2.5)


legend()
ylim([10^(-25) 10^(-22)])
xlim([n+1 N]);
ax=gca;
ax.FontSize=20;
ax.FontName='Times new roman';
xlabel('\Omega/2\pi   (Hz)');
ylabel('S(\Omega)   (Hz^{-1/2})');