%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXISWSH qdotdot = v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qdotdot] = system_SAE_III_qdotdot(q,u)

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

dq = [ q(3); q(4)];



%% STOIXEIA TIS DUNAMIKHS EXISWHS GIA TON ROMPOTIKO VRAXIONA ME 2 ARTHRWSEIS
%% PINAKAS H =  [ H11 H12 ; H21 H22]T
H11 =  m2*((Lc2^2)+(L1^2)) + (Lc1^2)*m1 + Iz2 + Iz1 + m2*2*(L1*Lc2*cos(q(2))); 
H12 = (Lc2^2)*m2 + Iz2 + L1*Lc2*m2*cos(q(2));
H21=H12;
H22 = (Lc2^2)*m2 + Iz2;

H = [ H11 H12 ; H21 H22];
H_1 = inv(H);

%% PINAKAS C = [ C11 C12 ; C21 C22]T
C11 =  -m2*L1*Lc2*sin(q(2))*q(4);
C12 = -m2*L1*Lc2*sin(q(2))*(q(3)+q(4));
C21 = m2*L1*Lc2*sin(q(2))*q(3);
C22 = 0;

C = [ C11 C12 ; C21 C22];

%% PINAKAS G = [g1 g2]T
g1 =  m2*Lc2*g*cos(q(1)+q(2)) + (m2*L1 + m1*Lc1)*g*cos(q(1));
g2 = m2*Lc2*g*cos(q(1)+q(2));

G = [g1; g2];

%% DIAFORIKH EXISWSH d^2q/dt^2 ME q = [q1 q2]T , dq/dt = [dq1/dt dq2/dt]T

qdotdot = H_1*(u - C*dq - G);
qdotdot = qdotdot';
end 
