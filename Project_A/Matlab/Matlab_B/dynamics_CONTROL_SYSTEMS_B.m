%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EPILUSH DIAFORIKHS EXISWSHS B
function dstate = dynamics_SAE_III_B(t,state,r)

t

%% INPUT FUNCTIONS
if (r(1) == 1)
    dr=0;
    d2r=0;
elseif(r(1)== 1.2)
    dr=1.2;
    d2r=0;
elseif(r(1)== 0.4)
    dr=0.4;
    d2r=0;
elseif(r(1)== 0.04)
    dr=0.04;
    d2r=0;
 end

%% N(s) FUNCTION
a=0.06;
eo=0.2;
e = state(1);

if(e<=-eo)
    N=1;
elseif(e>-eo && e<eo)
    N=a;
elseif(e>=eo)
    N=1;
end


%% DIFFERENCTIAL EQUATIONS OF STATES e =  [x1 x2]T
dstate(1)=state(2);
dstate(2)= d2r + dr  -state(2)  - 4*N*e;
dstate=dstate';

end
