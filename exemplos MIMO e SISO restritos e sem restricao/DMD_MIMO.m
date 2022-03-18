%% Funções de transferência entrada-saída

Ts = 1; %minutos

GG11 = tf([0.000298 0.002566],[1 -1.715 0.755],Ts);
GG11.outputdelay = 27;
g11 = step(GG11);
GG12 = tf([0.007106  0.007004],[1 -0.9714],Ts);
GG12.outputdelay = 23;
g12 = step(GG12);
GG13 = tf([-0.03856  0.009548 1.798e-5],[1 -0.9042 5.752e-7],Ts);
GG13.outputdelay = 15;
g13 = step(GG13);

GG21 = tf([-0.02409 0.0234],[1 -0.9984  0.02156],Ts);
g21 = step(GG21);
GG22 = tf([0],[1]);
GG23 = tf([0],[1]);

GG31 = tf([0.0009344],[1 -0.9445],Ts);
g31 = step(GG31);
GG32 = tf([-0.004995],[1 -0.9445],Ts);
g32 = step(GG32);
GG33 = tf([-0.0002795 0.0003478 0.0003198],[1 -1.652 0.6869],Ts);
GG33.outputdelay = 2;
g33 = step(GG33);

%% Funções de transferência perturbação

GQ11 = tf([0],[1]);
GQ12 = tf([0.01366 -0.0128 -0.01706],[1 -1.915 0.9191],Ts);
GQ12.outputdelay = 20;
q12 = step(GQ12);
GQ13 = tf([0],[1]);

GQ21 = tf([0.0004086  -0.0003754],[1 -1.306 0.3331],Ts);
q21 = step(GQ21);
GQ22 = tf([0.07873 -0.07549],[1 -1.176  0.2297],Ts);
q22 = step(GQ22);
GQ23 = tf([0],[1]);

GQ31 = tf([3.153e-5 -2.938e-5],[1 -1.9 0.903],Ts);
q31 = step(GQ31);
GQ32 = tf([0.008119],[1 -0.956],Ts);
q32 = step(GQ32);
GQ33 = tf([0],[1]);


%% Matriz G para simulação

GTFZ = [GG11 GG12 GG13 GQ11 GQ12 GQ13;
        GG21 GG22 GG23 GQ21 GQ22 GQ23;
        GG31 GG32 GG33 GQ31 GQ32 GQ33];



%% Calculo Nss

nss_io = zeros(3,3);
nss_q = zeros(3,2);
for i = 2:219
    if i <= 109
        if g11(i,1) ~= g11(i-1,1)
            nss_io(1,1) = i-1;
        end
    end
    if i <= 219
        if g12(i,1) ~= g12(i-1,1)
            nss_io(1,2) = i-1;
        end 
    end
    if i <= 115
        if g13(i,1) ~= g13(i-1,1)
            nss_io(1,3) = i-1;  
        end
    end
    if i <= 162
        if g21(i,1) ~= g21(i-1,1)
            nss_io(2,1) = i-1;
        end
    end
    if i <= 104
        if g31(i,1) ~= g31(i-1,1)
            nss_io(3,1) = i-1;  
        end
    end
    if i <= 104
        if g32(i,1) ~= g32(i-1,1)
            nss_io(3,2) = i-1; 
        end
    end    
    if i <= 47
        if g33(i,1) ~= g33(i-1,1)
            nss_io(3,3) = i-1;
        end
    end
    if i <= 162
        if q12(i,1) ~= q12(i-1,1)
            nss_q(1,2) = i-1;
        end
    end
    if i <= 162
        if q21(i,1) ~= q21(i-1,1)
            nss_q(2,1) = i-1;    
        end
    end
    if i <= 81
        if q22(i,1) ~= q22(i-1,1)
            nss_q(2,2) = i-1;    
        end  
    end
    if i <= 81
        if q31(i,1) ~= q31(i-1,1)
            nss_q(3,1) = i-1;    
        end 
    end
    if i <= 127
        if q32(i,1) ~= q32(i-1,1)
            nss_q(3,2) = i-1;    
        end
    end    
end

%% Definição dos parâmetros MPC
  
N1 = [15 1 1]';  
  
N2 = [50 52 47]';

N_1 = N2(1)-N1(1)+1;
N_2 = N2(2)-N1(2)+1;
N_3 = N2(3)-N1(3)+1;
sum_N = N_1+N_2+N_3;


Nu = [35  35  35]';

lambda_n = [1 1 1]';
delta_n = [1 1 1]';
lambda = [1 1 1]';
delta = [1 5 1]';

% lambda(1,1) = lambda_n(1,1)*(1^2*Nu(1,1))^-1;
% lambda(2,1) = lambda_n(2,1)*(1^2*Nu(2,1))^-1;
% lambda(3,1) = lambda_n(3,1)*(1^2*Nu(3,1))^-1;
% 
% delta(1,1) = delta_n(1,1)*(4^2*N_1)^-1;    %U1_max = 120
% delta(2,1) = delta_n(2,1)*(2^2*N_2)^-1;     %U2_max =23
% delta(3,1) = delta_n(3,1)*(2^2*N_3)^-1;    %U3_max = 190

I1 = eye(N2(1,1)-N1(1,1)+1);
I2 = eye(N2(2,1)-N1(2,1)+1);
I3 = eye(N2(3,1)-N1(3,1)+1);

Q_y = [delta(1,1)*I1 zeros(N2(1,1)-N1(1,1)+1,N2(2,1)-N1(2,1)+1) zeros(N2(1,1)-N1(1,1)+1,N2(3,1)-N1(3,1)+1);
       zeros(N2(2,1)-N1(2,1)+1,N2(1,1)-N1(1,1)+1) delta(2,1)*I2 zeros(N2(2,1)-N1(2,1)+1,N2(3,1)-N1(3,1)+1);
       zeros(N2(3,1)-N1(3,1)+1,N2(1,1)-N1(1,1)+1) zeros(N2(3,1)-N1(3,1)+1,N2(2,1)-N1(2,1)+1) delta(3,1)*I3]; 

I1 = eye(Nu(1,1));
I2 = eye(Nu(2,1));
I3 = eye(Nu(3,1));

Q_u = [lambda(1,1)*I1 zeros(Nu(1,1),Nu(2,1)) zeros(Nu(1,1),Nu(3,1));
       zeros(Nu(2,1),Nu(1,1)) lambda(2,1)*I2 zeros(Nu(2,1),Nu(3,1))
       zeros(Nu(3,1),Nu(1,1)) zeros(Nu(3,1),Nu(2,1)) lambda(3,1)*I3];
   


%% Matriz de ganhos G e H 

 G11 = [];
 G12 = [];
 G13 = [];
 G21 = [];
 G31 = [];
 G32 = [];
 G33 = [];
 aux = [];
 
for l = 1:Nu(1,1)
    aux = [zeros(l-1,1); g11(N1(1,1):(N2(1,1)-l+1))];
    G11 = [G11 aux];
    
    aux = [zeros(l-1,1); g21(N1(2,1):(N2(2,1)-l+1))];
    G21 = [G21 aux];
    
    aux = [zeros(l-1,1); g31(N1(3,1):(N2(3,1)-l+1))];
    G31 = [G31 aux];
end

for l = 1:Nu(2,1)
    aux = [zeros(l-1,1); g12(N1(1,1):(N2(1,1)-l+1))];
    G12 = [G12 aux];
   
    aux = [zeros(l-1,1); g32(N1(3,1):(N2(3,1)-l+1))];
    G32 = [G32 aux];
end

for l = 1:Nu(3,1)
    aux = [zeros(l-1,1); g13(N1(1,1):(N2(1,1)-l+1))];
    G13 = [G13 aux];
      
    aux = [zeros(l-1,1); g33(N1(3,1):(N2(3,1)-l+1))];
    G33 = [G33 aux];
end
G22 = zeros(N2(2,1)-N1(2,1)+1,Nu(2,1));
G23 = zeros(N2(2,1)-N1(2,1)+1,Nu(3,1));

G = [G11 G12 G13;
     G21 G22 G23;
     G31 G32 G33]; 
 
H11 = zeros(N2(1,1)-N1(1,1)+1,max(nss_io(:,1)));
H12 = zeros(N2(1,1)-N1(1,1)+1,max(nss_io(:,2)));
H13 = zeros(N2(1,1)-N1(1,1)+1,max(nss_io(:,3)));
H21 = zeros(N2(2,1)-N1(2,1)+1,max(nss_io(:,1)));
H22 = zeros(N2(2,1)-N1(2,1)+1,max(nss_io(:,2)));
H23 = zeros(N2(2,1)-N1(2,1)+1,max(nss_io(:,3)));
H31 = zeros(N2(3,1)-N1(3,1)+1,max(nss_io(:,1)));
H32 = zeros(N2(3,1)-N1(3,1)+1,max(nss_io(:,2)));
H33 = zeros(N2(3,1)-N1(3,1)+1,max(nss_io(:,3)));

for j = 1:max(nss_io(:,1))
    for i = N1(1,1):N2(1,1)
        if j <= size(g11,1) 
            if i+j <= nss_io(1,1)
                H11(i-N1(1,1)+1,j) = g11(i+j,1)-g11(j);
            else
                H11(i-N1(1,1)+1,j) = g11(nss_io(1,1),1)-g11(j);
            end    
        else
            H11(i-N1(1,1)+1,j) = 0;
        end
    end
    for i = N1(2,1):N2(2,1)
        if j <= size(g21,1) 
            if i+j <= nss_io(2,1)
                H21(i-N1(2,1)+1,j) = g21(i+j,1)-g21(j);
            else
                H21(i-N1(2,1)+1,j) = g21(nss_io(2,1),1)-g21(j);
            end    
        else
            H21(i-N1(2,1)+1,j) = 0;
        end
    end
    for i = N1(3,1):N2(3,1)
        if j <= size(g31,1) 
            if i+j <= nss_io(3,1)
                H31(i-N1(3,1)+1,j) = g31(i+j,1)-g31(j);
            else
                H31(i-N1(3,1)+1,j) = g31(nss_io(3,1),1)-g31(j);
            end    
        else
            H31(i-N1(3,1)+1,j) = 0;
        end
    end
end
 
for j = 1:max(nss_io(:,2))
    for i = N1(1,1):N2(1,1)
        if j <= size(g12,1) 
            if i+j <= nss_io(1,2)
                H12(i-N1(1,1)+1,j) = g12(i+j,1)-g12(j);
            else
                H12(i-N1(1,1)+1,j) = g12(nss_io(1,2),1)-g12(j);
            end    
        else
            H12(i-N1(1,1)+1,j) = 0;
        end

    end
    for i = N1(3,1):N2(3,1)
        if j <= size(g32,1) 
            if i+j <= nss_io(3,2)
                H32(i-N1(3,1)+1,j) = g32(i+j,1)-g32(j);
            else
                H32(i-N1(3,1)+1,j) = g32(nss_io(3,2),1)-g32(j);
            end    
        else
            H32(i-N1(3,1)+1,j) = 0;
        end
    end
end
for j = 1:max(nss_io(:,3))
    for i = N1(1,1):N2(1,1)
        if j <= size(g13,1) 
            if i+j <= nss_io(1,3)
                H13(i-N1(1,1)+1,j) = g13(i+j,1)-g13(j);
            else
                H13(i-N1(1,1)+1,j) = g13(nss_io(1,3),1)-g13(j);
            end    
        else
            H13(i-N1(1,1)+1,j) = 0;
        end
    end    
    for i = N1(3,1):N2(3,1)
        if j <= size(g33,1) 
            if i+j <= nss_io(3,3)
                H33(i-N1(3,1)+1,j) = g33(i+j,1)-g33(j);
            else
                H33(i-N1(3,1)+1,j) = g33(nss_io(3,3),1)-g33(j);
            end    
        else
            H33(i-N1(3,1)+1,j) = 0;
        end
    end    
end


H = [H11 H12 H13;
     H21 H22 H23;
     H31 H32 H33];

%% Calculo controlador

K = inv(G'*Q_y*G + Q_u)*G'*Q_y;
K1 = K(1,:);

%% Parametros simulação

t_sim = 1800;

%% inicialização variáveis

du_p = zeros(max(nss_io(:,1)) + max(nss_io(:,2)) + max(nss_io(:,3)),1);
du = zeros(Nu(1)+Nu(2)+Nu(3),1);
U = zeros(3,t_sim);
q = zeros(3,t_sim);
y = zeros(sum_N,1);
y_k = zeros(3,t_sim);
I = [ones(N2(1,1)-N1(1,1)+1,1) zeros(N2(1,1)-N1(1,1)+1,1) zeros(N2(1,1)-N1(1,1)+1,1);
     zeros(N2(2,1)-N1(2,1)+1,1) ones(N2(2,1)-N1(2,1)+1,1) zeros(N2(2,1)-N1(2,1)+1,1);
     zeros(N2(3,1)-N1(3,1)+1,1) zeros(N2(3,1)-N1(3,1)+1,1) ones(N2(3,1)-N1(3,1)+1,1)];
y_simulado = [];

t = 1:t_sim;
 
%% trajetorias de simulação;

r = [2*ones(1,t_sim);
     0*ones(1,200) 1*ones(1,1600);
     0*ones(1,600) 1*ones(1,1200)];
 
q = [zeros(1,1000) 150*ones(1,800);
     zeros(1,1300) 0.1*ones(1,500);
     zeros(1,1300) 0*ones(1,500)];
%% Simulaçãoo

for i = 2:t_sim
    
    % Calcular saída atual
%     [Y] = Calculo_Ym(y_k,U,i);
    inputs = [U' q'];
    y_simulado = lsim(GTFZ,inputs,t);
    y_k(:,i) = y_simulado(i,:)';   
    y = G*du + H*du_p + I*y_k(:,i);
    
    du_p = [du(1,1);
            du_p(1:160,1);
            du(Nu(1)+1,1);
            du_p(161+1:378,1);
            du(Nu(1)+Nu(2)+1,1);
            du_p(161+218+1:end-1,1)];    
        
    
    if i >= 1200
        r_a = [2*ones(1,N_1) ones(1,N_2) ones(1,N_3)]';
    else
        r_a = [r(1,i+N1(1,1):i+N2(1,1)) r(2,i+N1(2,1):i+N2(2,1)) r(3,i+N1(3,1):i+N2(3,1))]';
    end    
    
    du = K*(r_a-y);
    e = r_a-y;
    
    U(1,i) = U(1,i-1) + du(1,1);
    U(2,i) = U(2,i-1) + du(Nu(1)+1,1);
    U(3,i) = U(3,i-1) + du(Nu(1)+Nu(2)+1,1);
    
end   

%% plots
figure;
subplot(2,1,1);
plot(t,y_simulado(:,1))
hold on
plot(t,y_simulado(:,2))
plot(t,y_simulado(:,3))
plot(t,r(1,:))
plot(t,r(2,:))
plot(t,r(3,:))
legend('1','2','3')

hold off

subplot(2,1,2);
plot(t,U(1,:))
hold on
plot(t,U(2,:))
plot(t,U(3,:))
legend('1','2','3')



%%
U1 = U(:,1:end-1);
U2 = U(:,2:end);

delta_u = U2-U1;

plot(1:1799,delta_u);
legend('1','2','3');


%%

figure;
subplot(2,1,1);
plot(t,y_simulado(:,1))
hold on
plot(t,y_simulado(:,2))
plot(t,y_simulado(:,3))
plot(t,r(1,:))
plot(t,r(2,:))
plot(t,r(3,:))
legend('1','2','3')

hold off







