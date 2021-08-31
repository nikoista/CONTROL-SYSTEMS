%% SAE III EXETASTIKH IAN-FEB 2020-2021
%% NIKOLAOS ISTATIADIS  AEM:9175

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE PORTRAIT
function [] = phase_Portrait(dxdt,dydt)
[X,Y] = meshgrid(linspace(-4,4,40));
[sx,sy] = meshgrid(linspace(-2,2,10));

streamline(stream2(X, Y, ...                % Points
                   dxdt(X,Y), dydt(X,Y),... % Derivatives
                   sx, sy));                % Starting points
end