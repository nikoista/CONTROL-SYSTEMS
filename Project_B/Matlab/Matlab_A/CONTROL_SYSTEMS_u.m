%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXISWSH NOMOU ELEXOU u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u] = control_SAE_III_u(q,v)


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
dq = [q(3); q(4)];
v=v';

%% STOIXEIA TIS DUNAMIKHS EXISWHS GIA TON ROMPOTIKO VRAXIONA ME 2 ARTHRWSEIS
%% PINAKAS H =  [ H11 H12 ; H21 H22]T
H11 =  m2*((Lc2^2)+(L1^2)) + (Lc1^2)*m1 + Iz2 + Iz1 + m2*2*(L1*Lc2*cos(q(2))); 
H12 = (Lc2^2)*m2 + Iz2 + L1*Lc2*m2*cos(q(2));
H21=H12;
H22 = (Lc2^2)*m2 + Iz2;

H = [ H11 H12 ; H21 H22];

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

%% EXISOSH TOU NOMOU ELEGXOU APO u =  [u1 u2]T 

u = H*v + C*dq + G;

end

