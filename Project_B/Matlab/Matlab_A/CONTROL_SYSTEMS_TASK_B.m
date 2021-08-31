%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THEMA B.1 & C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAMIKOPOISH ME SKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOU SFALMATOS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELOPOIHSHS NA EXEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLOUS STO -10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KAI PROSWMOIWSEIS
nfig = 1;
x0=0;

%% THELOUME GRAMMIKOPOIHSH ME SKB POLOUS STO -10 OPOTE APO THN THEORITIKH 
%% ANALUSH GNWRIZOUME TA K1,K2 
K1 = 20;
K2 = 100;

%% ARXIKES TIMES q0 , qd0
q01  = -87;
q02 = 167;
dq01 = 0;
dq02 = 0;

%% ARXIKES SUNTHIKES GIA THN ODE 
state0 = [q01 , q02 , dq01 ,dq02 ];

%% XRONIKO DIASTHMA PROSOMOIWSEIS
tspan=0:0.001:10;

%% SUNARTHSH POU ULOPOIOUME TIS DIAFORIKES EXISWSEIS TON z1,z2,z3,z4
f=@(t,state)dynamics_SAE_III(t,state,K1,K2);
d=@(t,state)eventsfun_SAE_III(t,state,x0);
options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);
[t,state]=ode23s(f,tspan,state0,options);

q1 = state(:,1);
q2= state(:,2);

dq3 = state(:,3);
dq4 = state(:,4);


%% EPITHIMHTH TROXIA qd ME  qd = [q1d q2d]T KAI PARAGWGOI TOUS
i=1;
for tspan = t(1):0.001:t(end)
q1d(i)  = q1_desire(tspan); 
q2d(i)    = q2_desire(tspan);
dq1d(i)  = dq1_desire(tspan) ;
dq2d(i) = dq2_desire(tspan);
d2q1d(i)= d2q1_desire(tspan) ; 
d2q2d (i) = d2q2_desire(tspan);
i=i+1;
end

qd=[q1d ; q2d];
dqd=[dq1d ; dq2d];
d2qd=[d2q1d ; d2q2d];
qd = qd';
dqd = dqd';
d2qd = d2qd';


%% DHMIOURGOUME PINAKES ME TA q , dq/dt , d^2q/dt^2 
q = [q1 q2];
dq = [dq3 dq4];
v = d2qd -  K1*(dq - dqd) - K2*(q - qd);
d2q = v;

%% BRISKO XANA THN EISODO ELEGXOU
for i=1:1:length(q)
u(i,:) = control_SAE_III_u(state(i,:),v(i,:));
end

%% DHMIOURGOUME PINAKES ME TA  SFALMATA PARAKOLOUTHHSHS e, de/dt ,d^2e/dt^2 
e =  q - qd;
de = dq - dqd;
d2e = d2q - d2qd;

%% PLOTS GIA TIS APOKRISEIS TOU SUTHMATOS
figure(nfig)
plot(t,q1,t,q1d,t,e(:,1),'-.','LineWidth',1);
titleName = strcat("Position of Joint 1 figure",int2str(nfig));
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q1','q1(desire)' ,'Error_q1');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,q2,t,q2d,t,e(:,2),'-.','LineWidth',1);
titleName = strcat("Position of Joint 2 figure ",int2str(nfig));
title(titleName);
ylabel('$Degrees$','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q2','q2(desire)' ,'Error_q1');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,dq3,t,dq1d,t,de(:,1),'-.','LineWidth',1);
titleName = strcat("Speed of Joint 1 figure",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q1dot','q1dot(desire)' ,'Error_q1dot');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,dq4,t,dq2d,t,de(:,2),'-.','LineWidth',1);
titleName = strcat("Speed of Joint 2 figure",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q2dot','q2dot(desire)' ,'Error_q2dot');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,d2q(:,1) ,t,d2q1d,t,d2e(:,1),'-.','LineWidth',1);
titleName = strcat("Acceleration of Joint 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q1dot_dot','q1dot_dot(desire)' ,'Error_q1dot_dot');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,d2q(:,2) ,t,d2q2d,t,d2e(:,2),'-.','LineWidth',1);
titleName = strcat("Acceleration of Joint 2  figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('q2dot_dot','q2dot_dot(desire)' ,'Error_q2dot_dot');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,u(:,1),'LineWidth',1);
titleName = strcat(" Control Input 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('u1');
saveas(gcf,sprintf('%s.png', titleName));


nfig = nfig + 1;
figure(nfig)
plot(t,u(:,2 ),'LineWidth',1);
titleName = strcat(" Control Input 2 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('u2');
saveas(gcf,sprintf('%s.png', titleName));


nfig = nfig + 1;
figure(nfig)
plot(t,e,'LineWidth',1);
titleName = strcat(" Error  figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('Error');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,e(:,1),'red','LineWidth',1);
titleName = strcat(" Error 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('e1');
saveas(gcf,sprintf('%s.png', titleName));


nfig = nfig + 1;
figure(nfig)
plot(t,e(:,2),'red','LineWidth',1);
titleName = strcat(" Error 2 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('e2');
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig + 1;
figure(nfig)
plot(t,de,'LineWidth',1);
titleName = strcat(" Derivative of Error  figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('Derivative of Error');
saveas(gcf,sprintf('%s.png', titleName));


nfig = nfig + 1;
figure(nfig)
plot(t,de(:,1),'red','LineWidth',1);
titleName = strcat(" Derivative of Error 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('edot1');
saveas(gcf,sprintf('%s.png', titleName));


nfig = nfig + 1;
figure(nfig)
plot(t,de(:,2),'red','LineWidth',1);
titleName = strcat(" Derivative of Error 2 figure ",int2str(nfig));
title(titleName);
ylabel('$ $','Interpreter','latex','fontsize',20);
xlabel('$t$','Interpreter','latex','fontsize',12);
legend('edot2');
saveas(gcf,sprintf('%s.png', titleName));

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
