syms A1 Lf Rf kg Lg Rg Lm Rm kc Jt Bt kw s real

P = A1*kg/(Rf + Lf*s)
G = kc/((Lg+Lm)*s+Rg+Rm)
F = 1/(Jt*s+Bt)
H = kw

G1 = simplify(F*G*P/(1+F*G*H))
G2 = simplify(-F/(1+F*G*H))

%%
clear all
s = tf('s');

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;

P = A1*kg/(Rf + Lf*s)
G = kc/((Lg+Lm)*s+Rg+Rm)
F = 1/(Jt*s+Bt)
H = kw

G1 = minreal(F*G*P/(1+F*G*H))
G2 = minreal(-F/(1+F*G*H))

pole(G1)
zero(G1)
pole(G2)
zero(G2)

pzmap(G1,G2)

%%
A = [-Bt/Jt 0 kc/Jt ; 0 -Rf/Lf 0 ; -kw/(Lg+Lm) kg/(Lg+Lm) -(Rg+Rm)/(Lg+Lm)]
B = [0 -1/Jt ; A1/Lf 0 ; 0 0]
C = eye(3)
D = zeros(3,2)
eig(A)

sys1 = ss(A,B,C,D)
zpk(sys1)
%%

[A,B,C,D]=linmod('lab02')
sys2 = ss(A,B,C,D)
zpk(sys2)

%% 2.2
figure(1)
sim('lab02_2')
sim('lab02_2b')

subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2))
ylabel('\omega_{m1} (rad/s)')
title('Réponse du système M1 à un signal échelon de 100V')

subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3))
ylabel('i_{m1} (A)')

subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4))
ylabel('i_{f1} (A)')
xlabel('Temps (s)')

figure(2)
subplot(3,1,1)
plot(ScopeDataCp(:,1), ScopeDataCp(:,2))
ylabel('\omega_{m1} (rad/s)')
title('Réponse du système M1 à un signal échelon de 100V et une perturbation à 5s')

subplot(3,1,2)
plot(ScopeDataCp(:,1), ScopeDataCp(:,3))
ylabel('i_{m1} (A)')

subplot(3,1,3)
plot(ScopeDataCp(:,1), ScopeDataCp(:,4))
ylabel('i_{f1} (A)')
xlabel('Temps (s)')

%% 2.3 M2


A1 = 4;
Lf = 0;
Rf = 50;
kg = 5;
Lg = 0; 
Rg = 1;
Lm = 0;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;

sim('lab02_3')

figure(3)
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeDataL(:,1), ScopeDataL(:,2))
ylabel('\omega_{m} (rad/s)')
title('Réponse des systèmes M1 et M2 à un signal échelon de 100V')
legend('M1','M2', 'Location', 'Southeast')
legend('boxoff')


subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3), ScopeDataL(:,1), ScopeDataL(:,3))
ylabel('i_{m} (A)')
legend('M1','M2', 'Location', 'Northeast')
legend('boxoff')


subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4), ScopeDataL(:,1), ScopeDataL(:,4))
ylabel('i_{f} (A)')
xlabel('Temps (s)')
legend('M1','M2', 'Location', 'Southeast')
legend('boxoff')

%% 2.4

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0; 
kw = 0.5;

sim('lab02_4')

figure(4)
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeDataB(:,1), ScopeDataB(:,2))
ylabel('\omega_m (rad/s)')
title('Réponse des systèmes M1 et M3 à un signal échelon de 100V')
legend('M1','M3', 'Location', 'Southeast')
legend('boxoff')

subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3), ScopeDataB(:,1), ScopeDataB(:,3))
ylabel('i_m (A)')
legend('M1','M3', 'Location', 'Southeast')
legend('boxoff')

subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4), ScopeDataB(:,1), ScopeDataB(:,4))
ylabel('i_f (A)')
xlabel('Temps (s)')
legend('M1','M3', 'Location', 'Southeast')
legend('boxoff')

%% 3.1

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;
A2 = 20;
kt = 0.25;
kpre = 0.2423076923;

G3 = minreal(kpre*A2*G1/(1+A2*G1*kt));
%G4 = ???;
[A,B,C,D] = linmod('lab02_3_1');
sys = ss(A,B,C,D)
zpk(sys)

%% valeurs numeriques

A = [-Bt/Jt 0 kc/Jt ; -(A1*A2*kt)/Lf -Rf/Lf 0 ; -kw/(Lg+Lm) kg/(Lg+Lm) -(Rg+Rm)/(Lg+Lm)]
B = [0 -1/Jt ; kpre*A2*A1/Lf 0 ; 0 0]
C = eye(3)
D = zeros(3,2)
sys = ss(A,B,C,D)
zpk(sys)


%% 3.2 M5

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;
A2 = 20;
kt = 0.25;
kpre = 0.2423076923;

M1 = sim('lab02_2', 'StopTime', '3');
ScopeData = M1.get('ScopeData');
M4 = sim('lab02_3_2', 'StopTime', '3');
ScopeData3_2 = M4.get('ScopeData3_2');

figure(5)
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData3_2(:,1), ScopeData3_2(:,2))
ylabel('\omega_m (rad/s)')
title('Réponse des systèmes M1 et M4 à un signal échelon de 100V')
legend('M1','M4', 'Location', 'Southeast')
legend('boxoff')

subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3), ScopeData3_2(:,1), ScopeData3_2(:,3))
ylabel('i_m (A)')
legend('M1','M4', 'Location', 'Northeast')
legend('boxoff')

subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4), ScopeData3_2(:,1), ScopeData3_2(:,4))
ylabel('i_f (A)')
xlabel('Temps (s)')
legend('M1','M4', 'Location', 'Northeast')
legend('boxoff')

%% 3.2 (suite)

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;
A2 = 20;
kt = 0.25;
kpre = 0.2423076923;

M1 = sim('lab02_2b', 'StopTime', '10');
ScopeData = M1.get('ScopeDataCp');
M4 = sim('lab02_3_2b', 'StopTime', '10');
ScopeData3_2 = M4.get('ScopeData3_2');

figure(5)
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData3_2(:,1), ScopeData3_2(:,2))
ylabel('\omega_m (rad/s)')
title('Réponse des systèmes M1 et M4 à un signal échelon de 100V et perturbation à 5s')
legend('M1','M4', 'Location', 'Southeast')
legend('boxoff')

subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3), ScopeData3_2(:,1), ScopeData3_2(:,3))
ylabel('i_m (A)')
legend('M1','M4', 'Location', 'Northeast')
legend('boxoff')

subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4), ScopeData3_2(:,1), ScopeData3_2(:,4))
ylabel('i_f (A)')
xlabel('Temps (s)')
legend('M1','M4', 'Location', 'Northeast')
legend('boxoff')

%% 3.3

A1 = 4;
Lf = 2;
Rf = 50;
kg = 5;
Lg = 0.005; 
Rg = 1;
Lm = 0.005;
Rm = 1;
kc = 0.5;
Jt = 0.05;
Bt = 0.005; 
kw = 0.5;
A2 = 20;
kt = 0.25;
kpre = 0.2423076923;

M4 = sim('lab02_3_2', 'StopTime', '0.5');
ScopeData = M4.get('ScopeData3_2');

Lf = 0;
Lg = 0;
Lm = 0;

M5 = sim('lab02_3_2', 'StopTime', '0.5');
ScopeData3_2 = M5.get('ScopeData3_2');

figure(5)
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData3_2(:,1), ScopeData3_2(:,2))
ylabel('\omega_m (rad/s)')
title('Réponse transitoire des systèmes M4 et M5 à un échelon de 100V')
legend('M4','M5', 'Location', 'Southeast')
legend('boxoff')

subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,3), ScopeData3_2(:,1), ScopeData3_2(:,3))
ylabel('i_m (A)')
legend('M4','M5', 'Location', 'Northeast')
legend('boxoff')

subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,4), ScopeData3_2(:,1), ScopeData3_2(:,4))
ylabel('i_f (A)')
xlabel('Temps (s)')
legend('M4','M5', 'Location', 'Northeast')
legend('boxoff')
