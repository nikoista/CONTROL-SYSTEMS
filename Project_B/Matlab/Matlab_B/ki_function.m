%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUNARTHSH EURESHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ki DIANUSMATOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROSHMOY g(x)
function [k] = ki_function(q,qr)

n = [ 0.036 ; 0.036];

%% EKTIMOMENWN PARAMETROI ROMPOTIKOU VRAXIONA
L1 = 0.5;
L2 = 0.4;
Lc1 = 0.35;
Lc2 = 0.1;
m1 = 6;
m2 = 2;
Iz1 = 0.05;
Iz2 = 0.02;
g = 9.81;

%% PRAGMATIKOI PARAMETROI ROMPOTIKOU VRAXIONA
rL1 = 0.5;
rL2 = 0.4;
rLc1 = 0.2;
rLc2 = 0.4;
rm1 = 6;
rm2 = 4;
rIz1 = 0.43;
rIz2 = 0.05;
g = 9.81;


%% THETA 1 KAI THETA2 SUMFWNA ME TO max(|xekt - x| , |xekt - x|)
theta1 = [2.05; 3; 0.25;]; 
theta2 = [0.9925; 2.05; 3; 0.1125; 0.13; 0.45];


dq = [q(3) ; q(4)];
dqr = [qr(1) ; qr(2)];
d2qr = [qr(3) ; qr(4)];

w = q(4)*dqr(1) + (q(3) + q(4))*dqr(2);

Y1 = [-L1*sin(q(2))*w + g*cos(q(1) + q(2)) L1*g*cos(q(1)) m1*g*cos(q(1)); ...
    L1*sin(q(2))*dqr(1)*q(3) + g*cos(q(1) + q(2)) 0 0];

Y2 = [d2qr(1) + d2qr(2) , 2*L1*cos(q(2))*d2qr(1) +  L1*cos(q(2))*d2qr(2) , L1^2*d2qr(1) ,...
    m1*d2qr(1) , d2qr(1)+d2qr(2), d2qr(1);...
    d2qr(1) + d2qr(2) , L1*cos(q(2))*d2qr(1) , 0  , 0 , d2qr(1)+d2qr(2) , 0];

%% UPOLOGISMOS DIANUSMATOS k SUMFWNA ME THN THEORITIKH ANALYSH
k = abs(Y2)*theta2 + abs(Y1)*theta1 + n;


end

