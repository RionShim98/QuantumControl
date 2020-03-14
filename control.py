#例えば10^8を出したいときは10e7と打つ
#同様に10^(-2)を出したいときは10e-3と打つ

import matplotlib.pyplot as plt
import numpy as np
import math as ma
from numpy.linalg import inv
#from control import matlab

#調整可能なそれぞれのパラメータ
v=0#信号の分散
Qk=10*np.eye(10);Rk=0.001#Ricattiのパラメータ

#重力波検出器のパラメータ
hbar=1.055*10e-35;L=4000;M=40;Wm=1;P=800*10e2
ITM=2*ma.pi*200;c=3*10e7;lamlaser=1064*10e-10#レーザー周波数
Wo=2*ma.pi*c/lamlaser;Go=(Wo/L)*ma.sqrt(2*P*L/(hbar*Wo*c))
Gm=Go*ma.sqrt(hbar/(M*Wm))

#scaling law
Rsn=0.01;then=ma.pi/10
de=-ma.sqrt(Rsn)*ma.sin(2*then)*ITM/(1+Rsn+2*ma.sqrt(Rsn)*ma.cos(2*then))
IFO=(1-Rsn)*ITM/(1+Rsn+2*ma.sqrt(Rsn)*ma.cos(2*then))

#ロスのパラメータ
Tlos1=5*10^(-4);Tlos2=0;Tlos3=10^(-5);Tlos4=10^(-4)
L1=1.5;L2=L1;L3=300;Ln=0.5
los1=c*Tlos1/L1#a_1のロスレート
los2=c*Tlos2/L2#a_2のロスレート
los3=c*Tlos3/L3#a_3(MCC)のロスレート
los4=c*Tlos4/Ln#a_4(ループキャビティ)のロスレート

#MCCFBNDPAのパラメータ
kap1=2*c/L;lam=10^5;n1=2.01;gam=n1*lam
n2=0;kap2=n2*kap1;gNI=ma.sqrt(c*gam/(2*L))
g24=ma.sqrt(c*gam/Ln);g34=ma.sqrt(c*kap1/Ln)#Lnはループの長さ

#Original LIGOの状態行列
A1=np.matrix([[0,Wm,0,0],
            [0,0,ma.sqrt(2)*Gm,0],
            [0,0,-ITM/2,0],
            [ma.sqrt(2)*Gm,0,0,-ITM/2]])#[X P q p]
B1=np.matrix([[0,0,0],
            [ma.sqrt(1/(hbar*M*Wm)),0,0],
            [0,-ma.sqrt(ITM),0],
            [0,0,-ma.sqrt(ITM)]])#[F_{GW} Qin Pin]
C1=np.matrix([0,0,0,ma.sqrt(ITM)])#[X P q p]
D1=np.matrix([0,0,1])#[F_{GW} Qin Pin]

#最終周波数N,開始周波数n+1,
N=5*10**3;n=9;T=N-n#これによってn+1 Hz~N Hzのプロットが出来る
f=[];noiselist1=[];SQLlist=[]

for k in range(1,T):
    f.append(n+k)
    s=2*ma.pi*1j*(n+k)

    # Pout=G*[Qin Pin F_{GW}]
    G1=C1*inv(s*np.eye(4)-A1)*B1+D1
    ba=abs(G1[0,1]/G1[0,0]/(M*L*s**2))**2/2
    shot=abs(G1[0,2]/G1[0,0]/(M*L*s**2))**2/2
    noise1=ma.sqrt(ba+shot)
    noiselist1.append(noise1)

    SQL=ma.sqrt(abs(hbar/(M*L**2*s**2)))
    SQLlist.append(SQL)

plt.xscale("log");plt.yscale("log")
plt.plot(f,noiselist1)
plt.plot(f,SQLlist)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Sensitivity 1/rHz')
plt.ylim(10e-26,10e-23)
plt.xlim(n+1,N)
plt.show()

