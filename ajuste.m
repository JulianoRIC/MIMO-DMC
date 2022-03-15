function [N1, N2, N, Nu, t_5] = ajuste( FT, d, Ts)

h = stepinfo(FT, 'SettlingTimeThreshold', 0.05);
t_5 = h.SettlingTime;
N1 = d +1;      %hzonte pred mínimo (comeca apos atraso de d amostras)
N2 = ceil(median([((0.5*t_5)/Ts), (t_5/Ts)])); %hzonte pred max
N = (N2 - N1 + 1); %hzonte de predicao               
Nu= ceil(mean([N/5, N/2]));  %Horizonte de controle 

end

