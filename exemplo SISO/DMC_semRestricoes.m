clc
clear all
close all

%Algoritmo DMC SISO - Sem restrições

%% Modelo  do Sistema
Ts=0.5;                                %periodo de amostragem
G=tf([0 1],[0.1 1.1 1]);         %FT da planta continua
G.outputdelay=0.4;             %atraso
Gz =c2d(G,Ts);                   %FT da planta discretizada
d=1;                                      % atraso de tempo discreto
[B,A]=tfdata(Gz, 'v');   

ac = size(A,2); % nro de colunas de A
bc=  size(B,2);  %nro de colunas de B

%% Vector de coeficientes Gi
gi = step(Gz)';            %resposta ao degrau 

%gi =  [0 0.0355 0.3905 0.6301 0.7757 0.8639 0.9175 0.9499 0.9696 0.9816 ...
%       0.9888  0.9932  0.9959 0.9975 0.9985 0.9991 0.9994 0.9997 0.9998 ...
%       0.9999  0.9999  1.0000 1.0000 1.0000]'; 

Nss=length(gi)-1;   %Nss: horizonte do modelo

%% Define os parametros de sintonia do Controle Preditivo
t_5 = 3.57;     %tempo de assentamento
N1 = d +1;      %hzonte pred mínimo (comeca apos atraso de 1 amostra)
N2 = ceil(median([((0.5*t_5)/Ts), (t_5/Ts)])); %hzonte pred max
N = (N2 - N1 + 1); %hzonte de predicao               
Nu= ceil(mean([N/5, N/2]));  %Horizonte de controle 
lambda=1;           %ponderacao do controle
delta=1*eye(N2);   %Matriz de ponderacao do erro

%% Calculo da Matriz G
G = zeros(N2,Nu);
G(:,1) = gi(d+1:N2+d);
for i =2:Nu
   for j = 2:N2
    G(j,i) = G(j-1,i-1);
   end
end

%% calcula matriz K  e seleciono primeira linha
K =inv(G'*delta*G + lambda*eye(Nu))*G';
K1=  K(1,:);

%iteracoes
iters = 500;

%%%%%%%%%SIMULACAO%%%%%%%%%%%%%%%%
%% Condicoes iniciais
y_ant = zeros(1,ac); 
u_ant=zeros(1,bc);
u_livre = zeros(1,Nss);
u=zeros(1,iters);
deltaU=0;
 
%vetores para guardar saidas e perturbacao
 y = zeros(1,iters); %saida
 d=zeros(1,iters);  %perturbacao
 
 %referencia
 w = [zeros(1,2)  4*ones(1,200)  10*ones(1, iters) ]; 

 %% Loop control:
 
 for i=1:iters  
 
     y(i) = A*y_ant'+B*u_ant'; %saida
     
 %perturbacao degrau
 if i>400 
     d(i)=y(i)+1.5; 
 else 
     d(i) = y(i); 
 end
    
 %% Algoritmo DMC
 
 %resposta livre
 free=zeros(1,N2); 
 
    for(k=1:N2)
 
    for j=1:Nss-k 
        g_free(j) = gi(k+j)-gi(j);
    end
    for j=Nss-k+1:Nss
        g_free(j)= gi(Nss)-gi(j);
    end
    free(k) = d(i) + g_free*u_livre';
    end
    
    %considera ref futura cte
    ref=w(i)*ones(1,N2)';
    
    %calcula o controle 
    deltaU = K1*(ref-free'); %deltaU=k(w-f)
    if i == 1
        u(i)=deltaU;
    else
        u(i) = u(i-1) + deltaU;
    end
    
    %atualizacao do vetor de controle
   utemp = u_ant(1:bc-1);
   u_ant = [u(i)  utemp];
   
 end
 
%%  Graficos
subplot(3,1,1)
plot(d(1:iters), '--b');
hold 
plot(w(1:iters), 'm');
legend('saida', 'ref')

subplot(3,1,2)
p = [zeros(1,399) 0.5*ones(1,iters)];
plot(p(1:500), 'g');
legend('perturbacao')

subplot(3,1,3)
plot(u(1:iters),'-r');
legend('sinal de controle')



