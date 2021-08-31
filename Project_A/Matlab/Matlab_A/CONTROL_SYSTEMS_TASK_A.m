%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    THEMA A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N(s) DEN UPARXEI KAI   e=u


%% ARXIKOPOIHSH XRONIKOU DIASTHMATOS [0 ,20]
nfig = 0;
tspan=(0:0.1:20);

%% ARXIKOPOIHSH METABLHTWN KATASTASEIS

state0(1,:)=[-2 , 1.5];
state0(2,:)=[-2.5 , 0.8];
state0(3,:)=[1.5 , 2];
state0(4,:)=[0.2 , 1.8];
state0(5,:)=[2.5 , -0.8];
state0(6,:)=[2 , -2];
state0(7,:)=[-0.2 , -1.8];
state0(8,:)=[-1 , -2.5];
x0=0;

%% ARXIKOPOIHSH PINAKWN
x1_A = zeros(length(tspan),length(state0));
x2_A = zeros(length(tspan),length(state0));
error_A = zeros(length(tspan),length(state0));
x1dot_A = zeros(length(tspan),length(state0));
x2dot_A = zeros(length(tspan),length(state0));
l = zeros(length(tspan),length(state0));
x1_B = zeros(length(tspan),length(state0));
x2_B = zeros(length(tspan),length(state0));
error_B = zeros(length(tspan),length(state0));
x1dot_B = zeros(length(tspan),length(state0));
x2dot_B = zeros(length(tspan),length(state0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     a)

%% ARXIKOPOIHSH EISODOU r(t) = 1 BHMATIKH SUNARTHSH

r=@(t) 1;
dr=0;
d2r=0;

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
i

f=@(t,state)dynamics_SAE_III_A(t,state,r);
d=@(t,state)eventsfun_SAE_III(t,state,x0);

options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);

[t,state]=ode45(f,tspan,state0(i,:),options);


x1_A(:,i)=state(:,1);
x2_A(:,i)=state(:,2);

error_A(:,i) = state(:,1);

x1dot_A(:,i)=state(:,2);
x2dot_A(:,i)= d2r + dr -4*state(:,1)-state(:,2);


end

%% SUNARTHSH EISODOU - BHMATIKH SUNARTHSH
k=1;
for i=0:0.1:20
    l(k)=1;
    k=k+1;
end

%% PLOTS GIA TA x1,x2,e,x1dot,x2dot,r SE SUNARTHSH ME TON XRONO t
nfig = nfig +1;
figure(nfig);
plot(t,x1_A);
titleName = strcat("x1(t) - t ,  r(t)= 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ x1(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2_A);
titleName = strcat("x2(t) - t ,  r(t)= 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ x2(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,error_A);
titleName = strcat("e(t) - t ,  r(t)= 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ e(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x1dot_A);
titleName = strcat("x1dot - t , r(t)= 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ dx1(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2dot_A);
titleName = strcat("x2dot - t ,  r(t)= 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ dx2(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(tspan,l);
titleName = strcat("r(t) - t , r(t) = 1 figure ",int2str(nfig));
title(titleName);
ylabel('$ r(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
saveas(gcf,sprintf('%s.png', titleName));

%% PHASE PORTRAIT TWN x2,x1
size_of_states=size(x1_A);
nfig = nfig +1;
for i=1:size_of_states(2)
    figure(nfig);
    a =quiver(x1_A(:,i),x2_A(:,i),x1dot_A(:,i),x2dot_A(:,i))
    hold on;
    titleName = strcat("X1 - X2 ,  r(t)= 1 , N(s) = 0 (u=e) figure ",int2str(nfig));
    title(titleName);
    ylabel('$ X2 $','Interpreter','latex','fontsize',12)
    xlabel('$ X1 $','Interpreter','latex','fontsize',12)
end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));

%% PHASE PORTRAIT ME ALLH MORFH
nfig = nfig +1;
figure(nfig);
titleName = strcat("X1 - X2  figure",int2str(nfig));
title(titleName);
dxdt = @(x1,x2) x2;
dydt = @(x1,x2) -x2-4*x1+d2r+dr;

phase_Portrait(dxdt,dydt);
saveas(gcf,sprintf('%s.png', titleName));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     b)
%% ARXIKOPOIHSH EISODOU r(t) = 1.2*t SUNARTHSH RAMPA
r=@(t) 1.2*t;
dr=1.2;
d2r=0;

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
    
f=@(t,state)dynamics_SAE_III_A(t,state,r);
d=@(t,state)eventsfun_SAE_III(t,state,x0);

options=odeset('Events',d,'RelTol',10^-10,'AbsTol',10^-11);

[t,state]=ode45(f,tspan,state0(i,:),options);

x1_B(:,i)=state(:,1);
x2_B(:,i)=state(:,2);

error_B(:,i) = state(:,1);

x1dot_B(:,i)=state(:,2);
x2dot_B(:,i)= d2r +dr-4*state(:,1)-state(:,2);

end

%% PLOTS GIA TA x1,x2,e,x1dot,x2dot,r SE SUNARTHSH ME TON XRONO t
nfig = nfig +1;
figure(nfig);
plot(t,x1_B);
titleName = strcat("x1(t) - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ x1(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2_B);
titleName = strcat("x2(t) - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ x2(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,error_B);
titleName = strcat("e(t) - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ e(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x1dot_B);
titleName = strcat("x1dot - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx1(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2dot_B);
titleName = strcat("x2dot - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx2(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

%% SUNARTHSH EISODOU - RAMPA 1.2t
nfig = nfig +1;
figure(nfig);
plot(tspan,r(tspan));
titleName = strcat("r(t) - t , r(t)= 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ r(t) $','Interpreter','latex','fontsize',12);
xlabel('$ t $','Interpreter','latex','fontsize',12)
saveas(gcf,sprintf('%s.png', titleName));


%% PHASE PORTRAIT TWN x2,x1

size_of_states=size(x1_B);
nfig = nfig +1;
for i=1:size_of_states(2)
    figure(nfig);
    b = quiver(x1_B(:,i),x2_B(:,i),x1dot_B(:,i),x2dot_B(:,i))
    hold on;
    titleName = strcat("X1 - X2 ,  r(t)= 1.2t , N(s) = 0 (u=e) figure ",int2str(nfig));
    title(titleName);
    ylabel('$ X2 $','Interpreter','latex','fontsize',12)
    xlabel('$ X1 $','Interpreter','latex','fontsize',12)

end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));

%% PHASE PORTRAIT ME ALLH MORFH
nfig = nfig +1;
figure(nfig);
titleName = strcat("X1 - X2 figure",int2str(nfig));
title(titleName);
dxdt = @(x1,x2) x2;
dydt = @(x1,x2) -x2-4*x1+d2r+dr;

phase_Portrait(dxdt,dydt);
saveas(gcf,sprintf('%s.png', titleName));

%% SUGRHSH TWN 2 PHASE PORTRAIT

nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference N(s) = 0 (u=e) figure ",int2str(nfig));
title(titleName);
plot(x1_A,x2_A ,'red');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
plot(x1_B,x2_B,'blue');
legend({'Red Phase Port 1','Blue Phase Port 2'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));

%% PHASIKA PORTRAITA XWRIS QUIVER
nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference without N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_A,x2_A ,'red');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
legend({' Phase Port r=1'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));

%% PHASIKA PORTRAITA XWRIS QUIVER
nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference without N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_B,x2_B ,'blue');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
legend({' Phase Port r=1.2t'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));