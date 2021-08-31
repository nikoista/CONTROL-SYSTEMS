%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THEMA B.2 & C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SXEDIASH EUROSTOU NOMOU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ELEGXOU ME THN METHODO 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLISTHHSHS SFALMATOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KAI PROSWMOIWSEIS

nfig = 1;
x0=0;

%% PINAKAS L 
L = [10 0 ; 0 10];

%% ARXIKES TIMES q0 , qd0
q01  = -87;
q02 = 167;
dq01 = 0;
dq02 = 0;

%% ARXIKES SUNTHIKES GIA THN ODE 
state0 = [q01 , q02 , dq01 ,dq02];

%% XRONIKO DIASTHMA PROSOMOIWSEIS
tspan=0:0.001:6;


%% SUNARTHSH POU ULOPOIOUME TIS DIAFORIKES EXISWSEIS TON z1,z2,z3,z4
f=@(t,state)robust_dynamics_SAE_III(t,state,L);
d=@(t,state)eventsfun_SAE_III(t,state,x0);
options=odeset('Events',d,'RelTol',10^-4,'AbsTol',10^-5);
[t,state]=ode23s(f,tspan,state0);


%% DHMIOURGIA TWN METAVLHTWN KATASTASHS q , dq/dt
q1 = state(:,1);
q2 =  state(:,2);

dq3 = state(:,3);
dq4 = state(:,4);

q = [q1 q2];
dq = [dq3 dq4];
qall = [ q1 q2 dq3 dq4];


%% THA DHMIOURGHSOUME XANA TON NOMO ELEGXOU KAI GIA NA VROUME TA u , s , 
%% qdotqdot KAI NA UPOLOGISMOUE TO SFALMA PARAKOLOUTHHSHS KATHWS KAI NA 
%% FTIAXOUME TA APARAITHTA PLOTS

% %% ARXIKOPOIHSH PINAKWN 
% q1d = zeros(length(t),1);
% q2d = zeros(length(t),1);
% dq1d = zeros(length(t),1);
% dq2d = zeros(length(t),1);
% d2q1d = zeros(length(t),1);
% d2q2d = zeros(length(t),1);

%% EPITHIMHTH TROXIA qd ME  qd = [q1d q2d]T KAI PARAGWGOI TOUS
i=1;
for i= 1:1:length(t)
q1d(i)  = q1_desire(t(i)); 
q2d(i)    = q2_desire(t(i));
dq1d(i)  = dq1_desire(t(i)) ;
dq2d(i) = dq2_desire(t(i));
d2q1d(i)= d2q1_desire(t(i)) ; 
d2q2d (i) = d2q2_desire(t(i));
i=i+1;
end

qd=[q1d ; q2d];
dqd=[dq1d ; dq2d];
d2qd=[d2q1d ; d2q2d];
qd = qd';
dqd = dqd';
d2qd = d2qd';

%% DHMIOURGIA TOU DIANUSMATOS "EPITHUMHTHS TAXHTUTAS" dqr/dt
dqr  =  (dqd')  - L*((q - qd)');
d2qr =  (d2qd') - L*((dq - dqd)');
qrall = [dqr ; d2qr];
qrall = qrall';
dqr = dqr';


%% ARXIKOPOIHSH PINAKWN 
% ueq = zeros(length(q),2);
% k = zeros(length(q),2);
% s = zeros(length(q),2);
% g_x = zeros(length(q),2);
% u = zeros(length(q),2);
% d2q = zeros(length(q),2);

for i=1:1:length(q)
ueq(i,:) = robust_control_SAE_III_u(qall(i,:),qrall(i,:));
k(i,:) = ki_function(qall(i,:),qrall(i,:));
s(i,:)  = dq(i,:) -  dqr(i,:);
g_x = smooth_g_function(s(i,:),eps);
u(i,:) = ueq(i,:) - k(i,:).*(g_x');
U = u(i,:)';
d2q(i,:) = system_SAE_III_qdotdot(qall(i,:),U);
end

display(size(s))

%% DHMIOURGOUME PINAKES ME TA  SFALMATA PARAKOLOUTHHSHS e, de/dt ,d^2e/dt^2 
% a = qall - 
e =  q - qd;
de = dq - dqd;
d2e = d2q - d2qd;

%%%%%%%%%%%%%%% PLOTS GIA TOU ROMPOTIKOU VRAXIONA ME 2 BATHMOUS ELEUTHERIAS


%% PLOT GIA TIS EXISWSEIS KATASTASHS q TIS EPITHHMHTES TROXIES
%% KAI TIS PARAGWGOUS KATHOS KAI TO SFALMA PARAKOLOUTHHSHS MAZI

figure(nfig)
subplot(2,1,1);
plot(t,q1,t,q1d,t,e(:,1),'-.','LineWidth',1);
titleName = strcat("q1 - q1(desire) ||  error = q1- q1(desire)");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q1','q1(desire)' ,'q1(error)');

% nfig = nfig + 1;
% figure(nfig)
subplot(2,1,2);
plot(t,q2,t,q2d,t,e(:,2),'-.','LineWidth',1);
titleName = strcat("q2 - q2(desire) ||  error = q2- q2(desire)");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q2','q2(desire)' ,'q2(error)');


nfig = nfig + 1;
figure(nfig)
subplot(2,1,1);
plot(t,dq3,t,dq1d,t,de(:,1),'-.','LineWidth',1);
titleName = strcat("q1dot - q1dot(desire) || error = q1dot- q1dot(desire)");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('dq1/dt','dq1(desire)/dt' ,'dq1/dt(error)');

% nfig = nfig + 1;
% figure(nfig)
subplot(2,1,2);
plot(t,dq4,t,dq2d,t,de(:,2),'-.','LineWidth',1);
titleName = strcat("q2dot - q2dot(desire) || error = q2dot- q2dot(desire)");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('dq2/dt','dq2(desire)/dt' ,'dq2/dt(error)');


% nfig = nfig + 1;
% figure(nfig)
% plot(t,d2q(:,1) ,t,d2q1d,t,d2e(:,1),'-.','LineWidth',1);
% titleName = strcat("q1dot_dot - q1dot_dot(desire) ||  error = q1dot_dot - q1dot_dot(desire)");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('d^2q1/dt^2','d^2q1(desire)/dt^2' ,'d^2q1/dt^2(error)');
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(t,d2q(:,2) ,t,d2q2d,t,d2e(:,2),'-.','LineWidth',1);
% titleName = strcat("q2dot_dot - q2dot_dot(desire) ||  error = q2dot_dot - q2dot_dot(desire)");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('d^2q2/dt^2','d^2q2(desire)/dt^2' ,'d^2q2/dt^2(error)');
% 
% 

%% PLOTS GIA TON NOMO ELEGXOU u = [u1 , u2]T

nfig = nfig + 1;
figure(nfig)
subplot(2,1,1);
plot(t,u(:,1),'LineWidth',1);
titleName = strcat("u1 - t");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('u1');


% nfig = nfig + 1;
% figure(nfig)
subplot(2,1,2);
plot(t,u(:,2),'LineWidth',1);
titleName = strcat("u2 - t");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('u2');


%% PLOTS GIA TO "EPITHUMHTHS TAXHTUTAS" dqr/dt KAI d^2qr/d^2
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(dqr(:,1),'LineWidth',1);
% titleName = strcat("q1rdot - t ");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('q1rdot');
% 
% 
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(dqr(:,2),'LineWidth',1);
% titleName = strcat("q2rdot - t ");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('q2rdot');
% 


% nfig = nfig + 1;
% figure(nfig)
% plot(d2qr(:,1),'LineWidth',1);
% titleName = strcat("q1rdot_dot - t ");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('q1rdot_dot');
% 
% 
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(d2qr(:,2),'LineWidth',1);
% titleName = strcat("q2rdot_dot - t ");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('q2rdot_dot');


%% PLOTS GIA TA SFALMATA PARAKOLOUTHHSHS
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(t,e(:,1),'LineWidth',1);
% titleName = strcat("e1 - t");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('e1');
% 
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(t,e(:,2),'LineWidth',1);
% titleName = strcat("e2 - t");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('e2');
% 
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(t,de(:,1),'LineWidth',1);
% titleName = strcat("e1dot - t");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('e1dot');
% 
% 
% nfig = nfig + 1;
% figure(nfig)
% plot(t,de(:,2),'LineWidth',1);
% titleName = strcat("e2dot - t");
% title(titleName);
% ylabel('$Degrees$','Interpreter','latex','fontsize',20);
% xlabel('$t$','Interpreter','latex','fontsize',12);
% legend('e2dot');


%% PLOTS SFALMATA PARAKOLOUTHHSHS KAI EPIFANEIA OLISTHSHS

nfig = nfig + 1;
figure(nfig)
subplot(2,1,1);
plot(e(:,1),de(:,1),'LineWidth',1);
titleName = strcat("e - de/dt ");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('e1','de1/dt');

% nfig = nfig + 1;
% figure(nfig)
subplot(2,1,2);
plot(e(:,2),de(:,2),'LineWidth',1);
titleName = strcat("e - de/dt ");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('e2','de2/dt');



nfig = nfig + 1;
figure(nfig)
subplot(2,1,1);
plot(t,s(:,1),'LineWidth',1);
titleName = strcat("sliding Surface ");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('s1');
  
subplot(2,1,2);
plot(t,s(:,2),'LineWidth',1);
titleName = strcat("sliding Surface ");
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('s2');
ZL = zlim();
pl = line([-10 10] ,[ -1 1], [ZL(1), 0]);
pl.Color = 'red';
pl.LineStyle = '--';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUNARTHSH q1desire ,q1desire KAI OI PARAGWGOI TOUS 
function [q1d] = q1_desire(t)

if ( t <= 5)
    [q1d]= -90 + 50*(1-cos(0.63*t));
elseif( t > 5)
    [q1d] = 10;
end
end
function q2d = q2_desire(t)

if ( t <= 5)
    q2d = 170 - 60*(1-cos(0.63*t));
elseif( t > 5)
    q2d = 50;
end
end

%% 1 PARAGWGO qd
function dq1d = dq1_desire(t)

if ( t <= 5)
    dq1d =  31.5000*sin(0.63*t);
elseif( t > 5)
    dq1d = 0;
end
end
function dq2d = dq2_desire(t)

if ( t <= 5)
    dq2d = -37.8*sin(0.63*t);
elseif( t > 5)
    dq2d = 0;
end
end

%% 2 PARAGWGO qd
function d2q1d = d2q1_desire(t)

if ( t<= 5)
    d2q1d =  19.845*cos(0.63*t);
elseif( t > 5)
    d2q1d = 0;
end
end
function d2q2d = d2q2_desire(t)

if ( t <= 5)
    d2q2d = -23.814*cos(0.63*t);
elseif( t > 5)
    d2q2d = 0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
