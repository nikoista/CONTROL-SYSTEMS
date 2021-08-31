%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    THEMA B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N(s) UPARXEI 

%% ARXIKOPOIHSH XRONIKOU DIASTHMATOS [0 ,20]
tspan=(0:0.1:20);
nfig=0;

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
a=0.06;
eo=0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     a)

%% ARXIKOPOIHSH EISODOU r(t) = 1 BHMATIKH SUNARTHSH
r=@(t) 1;
dr=0;
d2r=0;

%% ARXIKOPOIHSH PINAKWN
x1_A = zeros(length(tspan),length(state0));
x2_A = zeros(length(tspan),length(state0));
error_A = zeros(length(tspan),length(state0));
x1dot_A = zeros(length(tspan),length(state0));
x2dot_A = zeros(length(tspan),length(state0));

x1_B = zeros(length(tspan),length(state0));
x2_B = zeros(length(tspan),length(state0));
error_B = zeros(length(tspan),length(state0));
x1dot_B = zeros(length(tspan),length(state0));
x2dot_B = zeros(length(tspan),length(state0));

x1_C = zeros(length(tspan),length(state0));
x2_C = zeros(length(tspan),length(state0));
error_C = zeros(length(tspan),length(state0));
x1dot_C = zeros(length(tspan),length(state0));
x2dot_C = zeros(length(tspan),length(state0));

x1_D = zeros(length(tspan),length(state0));
x2_D = zeros(length(tspan),length(state0));
error_D = zeros(length(tspan),length(state0));
x1dot_D = zeros(length(tspan),length(state0));
x2dot_D = zeros(length(tspan),length(state0));

l = zeros(length(tspan),length(state0));

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
    
    f=@(t,state)dynamics_SAE_III_B(t,state,r);
    d=@(t,state)eventsfun_SAE_III(t,state,x0);
    
    options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);
    
    [t,state]=ode45(f,tspan,state0(i,:),options);
    
    x1_A(:,i)=state(:,1);
    x2_A(:,i)=state(:,2);
    
    error_A(:,i) = state(:,1);
    
    e = error_A(:,i);
    size_of_states=size(e);
    for j=1:1:size_of_states(1)
        if(e(j)<=-eo)
            N=1;
        elseif(e(j)>-eo && e(j)<eo)
            N=a;
        elseif(e(j)>=eo)
            N=1;
        end
        NsA(j,i) = e(j)*N;
        x2dot_A(j,i) = d2r + dr -state(j,2)  -4*NsA(j,i);
    end
    x1dot_A(:,i) = state(:,2);

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
    hold on;
    quiver(x1_A(:,i),x2_A(:,i),x1dot_A(:,i),x2dot_A(:,i),1)
    titleName = strcat("X1 - X2 ,  r(t)= 1 , N(s) figure ",int2str(nfig));
    title(titleName);
    ylabel('$x1$','Interpreter','latex','fontsize',20)
    xlabel('$x1$','Interpreter','latex','fontsize',12)
end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     b)


%% ARXIKOPOIHSH EISODOU r(t) = 1.2*t RAMPA SUNARTHSH
r=@(t) 1.2*t;
dr=1.2;
d2r=0;
e=0;
state=0;

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
    f=@(t,state)dynamics_SAE_III_B(t,state,r);
    d=@(t,state)eventsfun_SAE_III(t,state,x0);
    
    options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);
    
    [t,state]=ode45(f,tspan,state0(i,:),options);
    
    x1_B(:,i)=state(:,1);
    x2_B(:,i)=state(:,2);
    
    error_B(:,i) = state(:,1);

    e = error_B(:,i);
    size_of_states=size(e);
    for(j=1:1:size_of_states(1))
        if(e(j)<=-eo)
            N=1;
        elseif(e(j)>-eo && e(j)<eo)
            N=a;
        elseif(e(j)>=eo)
            N=1;
        end
        NsB(j,i) = e(j)*N;
        x2dot_B(j,i) =d2r+ dr -4*NsB(j,i)-state(j,2);
    end
    x1dot_B(:,i) = state(:,2);
end


%% PLOTS GIA TA x1,x2,e,x1dot,x2dot,r SE SUNARTHSH ME TON XRONO t
nfig = nfig +1;
figure(nfig);
plot(t,x1_B);
titleName = strcat("x1(t) - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ x1(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2_B);
titleName = strcat("x2(t) - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ x2(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,error_B);
titleName = strcat("e(t) - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ e(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x1dot_B);
titleName = strcat("x1dot - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx1(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2dot_B);
titleName = strcat("x2dot - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx2(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(tspan,r(tspan));
titleName = strcat("r(t) - t , r(t) = 1.2t figure ",int2str(nfig));
title(titleName);
ylabel('$ r(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
saveas(gcf,sprintf('%s.png', titleName));

%% PHASE PORTRAIT TWN x2,x1
size_of_states=size(x1_B);
nfig = nfig +1;
for i=1:size_of_states(2)
    figure(nfig);
    hold on;
    quiver(x1_B(:,i),x2_B(:,i),x1dot_B(:,i),x2dot_B(:,i))
    titleName = strcat("X1 - X2 ,  r(t)= 1.2t , N(s) figure ",int2str(nfig));
    title(titleName);
    ylabel('$x1$','Interpreter','latex','fontsize',20)
    xlabel('$x2$','Interpreter','latex','fontsize',12)
end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     c)


%% ARXIKOPOIHSH EISODOU r(t) = 1 BHMATIKH SUNARTHSH
r=@(t) 0.4*t;
dr=0.4;
d2r=0;
e=0;
state=0;

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
    
    f=@(t,state)dynamics_SAE_III_B(t,state,r);
    d=@(t,state)eventsfun_SAE_III(t,state,x0);
    
    options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);
    
    [t,state]=ode45(f,tspan,state0(i,:),options);
    
    x1_C(:,i)=state(:,1);
    x2_C(:,i)=state(:,2);
    
    error_C(:,i) = state(:,1);
    
    e = error_C(:,i);
    size_of_states=size(e);
    for(j=1:1:size_of_states(1))
        if(e(j)<=-eo)
            N=1;
        elseif(e(j)>-eo && e(j)<eo)
            N=a;
        elseif(e(j)>=eo)
            N=1;
        end
        NsC(j,i) = e(j)*N;
        
        x2dot_C(j,i) = d2r+ dr -4*NsC(j,i) -state(j,2);
    end
    x1dot_C(:,i) = state(:,2);
    
end

%% PLOTS GIA TA x1,x2,e,x1dot,x2dot,r SE SUNARTHSH ME TON XRONO t
nfig = nfig +1;
figure(nfig);
plot(t,x1_C);
titleName = strcat("x1(t) - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ x1(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2_C);
titleName = strcat("x2(t) - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ x2(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,error_C);
titleName = strcat("e(t) - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ e(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x1dot_C);
titleName = strcat("x1dot - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx1(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2dot_C);
titleName = strcat("x2dot - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ dx2(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(tspan,r(tspan));
titleName = strcat("r(t) - t , r(t)= 0.4t figure ",int2str(nfig));
title(titleName);
ylabel('$ r(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
saveas(gcf,sprintf('%s.png', titleName));

%% PHASE PORTRAIT TWN x2,x1
size_of_states=size(x1_C);
nfig = nfig +1;
for i=1:size_of_states(2)
    figure(nfig);
    hold on;
    quiver(x1_C(:,i),x2_C(:,i),x1dot_C(:,i),x2dot_C(:,i))
    titleName = strcat("X1 - X2 ,  r(t)= 0.4t , N(s) figure ",int2str(nfig));
    title(titleName);
    ylabel('$x1$','Interpreter','latex','fontsize',20)
    xlabel('$x2$','Interpreter','latex','fontsize',12)
    
end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     d)


%% ARXIKOPOIHSH EISODOU r(t) = 1 BHMATIKH SUNARTHSH
r=@(t) 0.04*t;
dr=0.04;
d2r=0;
e=0;

%% LOOP GIA NA VROUME THN x1 KAI x2 KATHWS KAI TO e 
for i=1:1:length(state0)
    f=@(t,state)dynamics_SAE_III_B(t,state,r);
    d=@(t,state)eventsfun_SAE_III(t,state,x0);
    
    options=odeset('Events',d,'RelTol',10^-5,'AbsTol',10^-6);
    
    [t,state]=ode45(f,tspan,state0(i,:),options);
    
    
    x1_D(:,i)=state(:,1);
    x2_D(:,i)=state(:,2);
    
    error_D(:,i) = state(:,1);
    
    e = error_D(:,i);
    size_of_states=size(e);
    for j=1:1:size_of_states(1)
        if(e(j)<=-eo)
            N=1;
        elseif(e(j)>-eo && e(j)<eo)
            N=a;
        elseif(e(j)>=eo)
            N=1;
        end
        NsD(j,i) = e(j)*N;
        x2dot_D(j,i) = d2r+ dr -4*NsD(j,i)-state(j,2);
    end
    x1dot_D(:,i) = state(:,2);
end

%% PLOTS GIA TA x1,x2,e,x1dot,x2dot,r SE SUNARTHSH ME TON XRONO t
nfig = nfig +1;
figure(nfig);
plot(t,x1_A);
titleName = strcat("x1(t) - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ x1(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2_A);
titleName = strcat("x2(t) - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ x2(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,error_A);
titleName = strcat("e(t) - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ e(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x1dot_A);
titleName = strcat("x1dot - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ dx1(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(t,x2dot_A);
titleName = strcat("x2dot - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ dx2(t)/dt $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
saveas(gcf,sprintf('%s.png', titleName));

nfig = nfig +1;
figure(nfig);
plot(tspan,r(tspan));
titleName = strcat("r(t) - t , r(t) = 0.04t  figure ",int2str(nfig));
title(titleName);
ylabel('$ r(t) $','Interpreter','latex','fontsize',12)
xlabel('$ t $','Interpreter','latex','fontsize',12)
saveas(gcf,sprintf('%s.png', titleName));


%% PHASE PORTRAIT TWN x2,x1
size_of_states=size(x1_D);
nfig = nfig +1;
for i=1:size_of_states(2)
    figure(nfig);
    hold on;
    quiver(x1_D(:,i),x2_D(:,i),x1dot_D(:,i),x2dot_D(:,i))
    titleName = strcat("X1 - X2 ,  r(t) = 0.04t , N(s) figure ",int2str(nfig));
    title(titleName);
    ylabel('$x1$','Interpreter','latex','fontsize',20)
    xlabel('$x2$','Interpreter','latex','fontsize',12)

end
legend({'(-2,1.5)','(-2.5,0.8)','(1.5,2)' , '(0.2,1.8)' , '(2.5,-0.8)','(2,-2)','(-0.2,1.8)','(-1,2.5)'})
hold off;

%% SUGRHSH TWN 2 PHASE PORTRAIT

nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference with N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_A,x2_A ,'red');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
plot(x1_A(1),x2_A(1),'s','LineWidth',2)
plot(x1_A(end),x2_A(end),'x','LineWidth',2)
legend({'Red Phase Port r=1'})%'Blue Phase Port r=1.2t','Blue Phase Port r=0.4t','Blue Phase Port r=0.04t'})
hold off;
saveas(gcf,sprintf('%s.png', titleName));


%%
nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference with N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_B,x2_B ,'blue');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
legend({' Phase Port r=1.2t'})
plot(x1_B(1),x2_B(1),'s','LineWidth',2)
plot(x1_B(end),x2_B(end),'x','LineWidth',2)
hold off;
saveas(gcf,sprintf('%s.png', titleName));


%%
nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference with N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_C,x2_C ,'green');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
legend({'Phase Port r=0.4t'})
plot(x1_C(1),x2_C(1),'s','LineWidth',2)
plot(x1_C(end),x2_C(end),'x','LineWidth',2)
hold off;
saveas(gcf,sprintf('%s.png', titleName));


%%
nfig = nfig +1;
figure(nfig);
hold on;
titleName = strcat("Phase Portraits Difference with N(s) figure ",int2str(nfig));
title(titleName);
plot(x1_D,x2_D ,'black');
ylabel('$ X2 $','Interpreter','latex','fontsize',12)
xlabel('$ X1 $','Interpreter','latex','fontsize',12)
legend({'Phase Port r=0.04t'})
plot(x1_D(1),x2_D(1),'s','LineWidth',2)
plot(x1_D(end),x2_D(end),'x','LineWidth',2)
hold off;
saveas(gcf,sprintf('%s.png', titleName));
