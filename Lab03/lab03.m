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

L = 0.00025;
R = 0.5;
Kc = 0.48;
Km = 0.5;
Bm = 0.02;
Jm = 0.05;

A=[-R/L -Kc/L ; Km/Jm -Bm/Jm];
B=[1/L ; 0];
C=[0 1];
D=[0];
I=eye(2);

F = C*inv(s*I-A)*B+D;
fprintf('ZPK:\n:')
zpk(F)
fprintf('Pole:\n:')
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
syms thetae B1 J1 l1 s

L = 0.00025
R = 0.5
Kc = 0.48
Km = 0.5
Bm = 0.02
Jm = 0.05
%l1 = 0.3
M1 = 2
%B1 = 0.125
g = 9.81
%J1 = 0

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
%% 3.1
syms Bm Jm
R = 0.5;
Km = 0.5;
Kc = 0.48;
%Fa = Km/(Bm*R + Kc*Km + Jm*R*s)
tau = interp1(ScopeData(:,2), ScopeData(:,1), 0.63*2)

plot(ScopeData(:,1), ScopeData(:,2), tau, 0.63*2, 'b*')
title('Gain statique et temps de montée')
xlabel('Temps (s)')
ylabel('\omega (rad/s)')

%% 3.2
delta_theta_max = max(theta(:,2))
Tp = theta((find(theta(:,2) == delta_theta_max)),1)
delta_theta_rp = theta(end)
delta_u = 0.017298519974968

plot(theta(:,1), theta(:,2))
title('Réponse du système en équilibre à une perturbation de 1% en entrée.')
xlabel('Temps (s)')
ylabel('\omega (rad/s)')

%% 4.1 Schema de simulation

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
J1 = 0.06
thetae = pi/4

A=[-R/L 0 -Kc/L
    0 0 1
    Km/(Jm+J1) -M1*g*l1*thetae*cos(thetae)/(2*(Jm+J1)) -(Bm+B1)/(Jm+J1)]
B = [1/L ; 0 ; 0]
C = [0 1 0]
D = [0]

%% 4.2 Perturbation condition initiale
perturbation_u = 0;

figure(1)
perturbation_theta = -30;
sim('Schema_de_simulation_42')
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_theta = -15;
sim('Schema_de_simulation_42')
subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_theta = -5;
sim('Schema_de_simulation_42')
subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

figure(2)
perturbation_theta = 5;
sim('Schema_de_simulation_42')
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_theta = 15;
sim('Schema_de_simulation_42')
subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_theta = 30;
sim('Schema_de_simulation_42')
subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
titre = sprintf('Réponse du système à une perturbation de %d degrés par rapport au point d''équilibre.', perturbation_theta)
title(titre)
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

thetae = pi/4;
A=[-R/L 0 -Kc/L
    0 0 1
    Km/(Jm+J1) -M1*g*l1*thetae*cos(thetae)/(2*(Jm+J1)) -(Bm+B1)/(Jm+J1)]
B = [1/L ; 0 ; 0]
C = [0 1 0]
D = [0]

%% 4.3 Perturbation entree
perturbation_theta = 0;
perturbation_u = 0.05;

figure(1)
perturbation_u = 0.05;
sim('Schema_de_simulation_42')
subplot(3,1,1)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
title('Réponse du système à une perturbation de 5% par rapport à la tension d''équilibre.')
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_u = 0.1;
sim('Schema_de_simulation_42')
subplot(3,1,2)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
title('Réponse du système à une perturbation de 10% par rapport à la tension d''équilibre.')
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

perturbation_u = 0.25;
sim('Schema_de_simulation_42')
subplot(3,1,3)
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
title('Réponse du système à une perturbation de 25% par rapport à la tension d''équilibre.')
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');

figure(2)
perturbation_u = 0.45;
sim('Schema_de_simulation_42')
plot(ScopeData(:,1), ScopeData(:,2), ScopeData(:,1), ScopeData(:,3))
title('Réponse du système à une perturbation de 45% par rapport à la tension d''équilibre.')
legend('Non-linéaire', 'Linéarisé', 'Location', 'best');