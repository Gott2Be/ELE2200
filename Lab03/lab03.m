syms L Kc R Km Bm Jm s

A=[-R/L -Kc/L ; Km/Jm -Bm/Jm]
B=[1/L ; 0]
C=[0 1]
D=[0]
I=eye(2)

F = simplify(C*inv(s*I-A)*B+D)

%%
clear all
s = tf('s')

L = 0.00025
R = 0.5
Kc = 0.48
Km = 0.5
Bm = 0.02
Jm = 0.05

A=[-R/L -Kc/L ; Km/Jm -Bm/Jm]
B=[1/L ; 0]
C=[0 1]
D=[0]
I=eye(2)

F = C*inv(s*I-A)*B+D
zpk(F)
pole(F)

%%

syms Kc R Km Bm Jm s

Fa = Km/(Bm*R + Kc*Km + Jm*R*s)

%%
clear all
s = tf('s')

L = 0.00025
R = 0.5
Kc = 0.48
Km = 0.5
Bm = 0.02
Jm = 0.05

A=[-R/L -Kc/L ; Km/Jm -Bm/Jm]
B=[1/L ; 0]
C=[0 1]
D=[0]
I=eye(2)

F = C*inv(s*I-A)*B+D
Fa = Km/(Bm*R + Kc*Km + Jm*R*s)

zpk(F)
pole(F)
zpk(Fa)
pole(Fa)

%% 2.2 moteur + bras
clear all
%s = tf('s')
syms thetae s

L = 0.00025
R = 0.5
Kc = 0.48
Km = 0.5
Bm = 0.02
Jm = 0.05
l1 = 0.3
M1 = 2
B1 = 0.125
g = 9.81
J1 = 0

A=[-R/L 0 -Kc/L
    0 0 1
    Km/(Jm+J1) -M1*g*l1*pi*cos(thetae)/(10*(Jm+J1)) -(Bm+B1)/(Jm+J1)]
B = [1/L ; 0 ; 0]
C = [0 1 0]
D = [0]
I=eye(3)

G =  simplify(C*inv(s*I-A)*B+D)

L = 0
A=[-R 0 -Kc
    0 0 1
    Km/(Jm+J1) -M1*g*l1*pi*cos(thetae)/(10*(Jm+J1)) -(Bm+B1)/(Jm+J1)]
B = [1 ; 0 ; 0]
C = [0 1 0]
D = [0]
I=eye(3)

Ga = simplify(C*inv(s*I-A)*B+D)
%% 3
syms Bm Jm
R = 0.5;
Km = 0.5;
Kc = 0.48;
Fa = Km/(Bm*R + Kc*Km + Jm*R*s)
interp1(ScopeData(:,2), ScopeData(:,1), 0.63*2)

%% 




