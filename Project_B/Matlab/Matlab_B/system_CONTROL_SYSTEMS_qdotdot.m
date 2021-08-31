%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXISWSH qdotdot = v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qdotdot] = system_SAE_III_qdotdot(q,u)

dq = [ q(3); q(4)];

%% PRAGMATIKOI PARAMETROI ROMPOTIKOU VRAXIONA
L1 = 0.5;
L2 = 0.4;
Lc1 = 0.2;
Lc2 = 0.4;
m1 = 6;
m2 = 4;
Iz1 = 0.43;
Iz2 = 0.05;
g = 9.81;

%% SLOTINE Example 4.9 Two link manipulator
theta1 = Iz1 + m1*(Lc1^2) + Iz2 + m2*(Lc2^2) + m2*(L1^2);
theta2 = Iz2 + m2*(Lc2^2);
theta3 = m2*L1*Lc2;
theta4 = m2*Lc2*g;
theta5 = (m2*L1 + m1*Lc1)*g;

H = [theta1+2*theta3*cos(q(2))  ,  theta2+theta3*cos(q(2)) ; theta2+theta3*cos(q(2)) , theta2];
C = [-theta3*sin(q(2))*(q(4))  ,  -theta3*sin(q(2))*(q(3) + q(4)) ; theta3*sin(q(2))*(q(3)) , 0];
G = [ theta4*cos(q(1) + q(2)) + theta5*cos(q(1)) ; theta4*cos(q(1) + q(2))];

H_1 = inv(H);

%% DIAFORIKH EXISWSH d^2q/dt^2 ME q = [q1 q2]T , dq/dt = [dq1/dt dq2/dt]T

qdotdot = H_1*(u - C*dq - G);
qdotdot = qdotdot';
end 
