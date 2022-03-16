function [ ss ] = steady(FT)

%steady-state values dos ensaios de malha aberta
ss = stepinfo(FT, 'SettlingTimeThreshold', 0.02);
ss = ss.SettlingMax;

end

