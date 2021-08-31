%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMALH SUNARTHSH 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROSHMOY g(x)
function [g] = smooth_g_function(x,eps)

x1= x(1);
x2= x(2);

if( abs(x1) >= eps)
    g(1) = x1/abs(x1);
else
    g(1)= x1/eps;
end

if( abs(x2) >= eps)
    g(2) = x2/abs(x2);
else
    g(2)= x2/eps;
end

g=g';
end