%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXISWSH EUROSTOU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOMOU ELEXOU u
function [ueq] = robust_control_SAE_III_u(q,qr)

dq = [q(3); q(4)];
dqr = [qr(1) ; qr(2)];
d2qr = [qr(3) ; qr(4)];

%% EKTIMOMENOI PARAMETROI ROMPOTIKOU VRAXIONA
L1 = 0.5;
L2 = 0.4;
Lc1 = 0.35;
Lc2 = 0.1;
m1 = 6;
m2 = 2;
Iz1 = 0.05;
Iz2 = 0.02;
g = 9.81;

%% DHMIOURGIA PINAKA GIA PIO EUKOLES PRAXEIS

theta1 = Iz1 + m1*(Lc1^2) + Iz2 + m2*(Lc2^2) + m2*(L1^2);
theta2 = Iz2 + m2*(Lc2^2);
theta3 = m2*L1*Lc2;
theta4 = m2*Lc2*g;
theta5 = (m2*L1 + m1*Lc1)*g;

H = [theta1+2*theta3*cos(q(2))  ,  theta2+theta3*cos(q(2)) ; theta2+theta3*cos(q(2)) , theta2];
C = [-theta3*sin(q(2))*(q(4))  ,  -theta3*sin(q(2))*(q(3) + q(4)) ; theta3*sin(q(2))*(q(3)) , 0];
G = [ theta4*cos(q(1) + q(2)) + theta5*cos(q(1)) ; theta4*cos(q(1) + q(2))];

%% EXISOSH TOU NOMOU ELEGXOU APO u =  [u1 u2]T 

ueq = H*d2qr + C*dqr + G;

end
