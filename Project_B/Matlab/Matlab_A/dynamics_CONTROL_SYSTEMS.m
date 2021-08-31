%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIAFORIKH EXISWSH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KATASTASEWN z1,z2,z3,z4
function [dz] = dynamics_SAE_III(t,q,K1,K2)

t

%% METAVLHTES KATASTASEIS q = [q1 q2]T  KAI PARAGWGOI TOUS
z(1) = q(1); %q1
z(2) = q(2); %q2
z(3) = q(3); %q1dot
z(4) = q(4); %q2dot

dz(1) = z(3);
dz(2) = z(4);


%% DHMIOURGIA PINAKWN GIA PIO EUKOLES PRAXEIS
nq = [q(1) q(2)];
qd = [q1_desire(t) q2_desire(t)];

dq = [q(3) q(4)];
dqd = [dq1_desire(t) dq2_desire(t)];
d2qd = [d2q1_desire(t) d2q2_desire(t)];

%% GRAMMIKOPOIHSH ME ANADRASH OPOU v = d^2q/dt^2
v = d2qd -  K1*(dq - dqd) - K2*(nq - qd);

%% NOMOS ELEXOU POU ANALYSAME STO ERWTHMA A
u = control_SAE_III_u(q,v);

%% UPOLOGISMOS TOU d^2q/dt^2 MESW TOU NOMOU ELEGXOU POU ANAPTUXAME
qdotdot = system_SAE_III_qdotdot(q,u);

dz(3) =  qdotdot(1);
dz(4) =  qdotdot(2);
dz=dz';

end


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