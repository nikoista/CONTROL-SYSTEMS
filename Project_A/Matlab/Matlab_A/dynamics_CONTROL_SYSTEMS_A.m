%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EPILUSH DIAFORIKHS EXISWSHS A
function dstate = dynamics_SAE_III_A(t,state,r)

t

%% INPUT FUNCTIONS
if (r(t) == 1)
    dr=0;
    d2r=0;
else
    dr=1.2;
    d2r=0;
end

%% DIFFERENCTIAL EQUATIONS OF STATES [x1 x2]T
dstate(1)=state(2);
dstate(2)= d2r+dr -state(2)-4*state(1);
dstate=dstate';

end

