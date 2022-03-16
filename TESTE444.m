clc
clear all
close all

%Algoritmo DMC MIMO - Sem restrições

%% Modelo  do Sistema
TS = 1;                                %periodo de amostragem
z    =  tf('z', TS);

%funcoes de transferencia em relacao a MV
G11 =  tf([0.000298 0.002566], [1 -1.751 0.755], TS, 'InputDelay', 27);
G12 =  tf([0.007106 0.007004], [1 -0.9714], TS, 'InputDelay', 23); 
G18 =  tf([-0.03856 0.009548 1.798e-5], [1 -0.9042 5.752e-7], TS, 'InputDelay', 15);
G25 =  tf([-0.022697 -0.0136], [1 -1.112 0.2116], TS); 
G35 =  tf([0.007741 0.003937], [1 -0.9584 1.897e-5], TS, 'InputDelay', 6);
G38 =  tf([0 -0.007009], [1 -0.9584], TS, 'InputDelay', 3);
G41 =  tf([-0.04949 0.01427 3.437e-5], [1 -0.968 1.17e-6], TS, 'InputDelay', 28);
G42 =  tf([-0.03625 0.007826], [1 -0.9672], TS, 'InputDelay', 23);
G44 =  tf([0 5.352], [1  -0.5134], TS);
G48 =  tf([0.02216 -0.01418], [1 -1.679 0.6918], TS, 'InputDelay', 14);
G54 =  tf([-0.2807 0.01537], [1 -1.469 0.5332], TS, 'InputDelay', 12); 
G56 =  tf([-0.04508 0.03108], [1 -1.249 0.3724], TS, 'InputDelay', 9);
G57 =  tf([0.04439 -0.03667], [1 -1.817 0.8272], TS, 'InputDelay', 5);
G66 =  tf([0.2632 0.1493], [1 -0.7659  0.1859], TS);
G71 =  tf([-0.02409 0.0234], [1 -0.9984 0.02156], TS);
G81 =  tf([0 0.0009344], [1 -0.9445], TS);
G82 =  tf([0 -0.004995], [1 -0.9445], TS); 
G83 =  tf([-0.006484 0.004228], [1 -1.432 0.4647], TS);
G88 =  tf([-0.0002795 0.0003478 0.0003198], [1 -1.652 0.6869], TS, 'InputDelay', 2);

%valores de regime permanente
ss_g11 = steady(G11);
ss_g12 = steady(G12);
ss_g18 = steady(G18);
ss_g25 = steady(G25);
ss_g35 = steady(G35);
ss_g38 = steady(G38);
ss_g41 = steady(G41);
ss_g42 = steady(G42);
ss_g44 = steady(G44);
ss_g48 = steady(G48);
ss_g54 = steady(G54);
ss_g56 = steady(G56);
ss_g57 = steady(G57);
ss_g66 = steady(G66);
ss_g71 = steady(G71);
ss_g81 = steady(G81);
ss_g82 = steady(G82);
ss_g83 = steady(G83);
ss_g88 = steady(G88);

ss = [ss_g11 ss_g12 ss_g18 ss_g25 ss_g35 ss_g38 ...
         ss_g41 ss_g42 ss_g44 ss_g48 ss_g54   ss_g56 ...  
         ss_g57 ss_g66 ss_g71 ss_g81  ss_g82 ss_g83 ss_g88]; 

%Matriz de funcoes transferencia em relacao as MV
Gu =  [
          G11 G12 0 0 0 0    0    G18; 
          0  0    0   0    G25   0   0   0;
          0  0  0  0    G35   0  0   G38;
          G41 G42 0 G44  0 0 0 G48;
          0 0 0 G54 0    G56  G57   0;
          0    0   0   0   0    G66   0   0;
          G71    0   0   0  0    0    0   0;
          G81 G82 G83 0  0 0  0  G88
          ];

%funcoes de transferencia em relacao a perturbacao
Gd12 =  tf([0.01366 -0.0128 -0.01706], [1 -1.915 0.9191], TS, 'InputDelay', 20);
Gd21   =  tf([0.00236 0.002173], [1 -0.9518 0.1157], TS);
Gd22   =  tf([4.68 -3.722], [1 -0.8786 3.875e-5], TS, 'InputDelay', 20);
Gd42 =  tf([-0.00611 0.005808 0.007654], [1 -1.914 0.9181], TS, 'InputDelay', 20);
Gd53   =  tf([0.03646 -0.006099 -0.02366], [1 -1.444 0.4712], TS, 'InputDelay', 6);
Gd63   =  tf([0.5015 0.2255], [1 -0.3526 0.0961], TS);
Gd71 =  tf([0.0004086 -0.0003754], [1 -1.306 0.3331], TS);
Gd72   =  tf([0.07873 -0.07549], [1 -1.176 0.2297], TS);
Gd81 =  tf([3.153e-5 -2.93810e-5], [1 -1.9  0.903], TS);
Gd82 =  tf([0 0.008119], [1 -0.956], TS);

% matriz perturbacao
Gw = [
           0         Gd12     0;
           Gd21  Gd22    0;
           0         0            0;
           0         Gd42     0;
           0         0            Gd53;
           0         0            Gd63;
           Gd71 Gd72     0;
           Gd81 Gd82     0
           ];
                
%% Vector de coeficientes Gi
%Nss: horizonte do modelo
g11  =  step(G11)'; 
g11(1) = [ ];
Nss_1 =  307;   
g12     = step(G12)';
Nss_2 = 197;  
g18     = step(G18)';
Nss_3 =  63;  
g25     = step(G25)';
g25(1) = [ ];
Nss_4 = 33;  
g35     = step(G35)';
g35(1) = [ ];
Nss_5 = 85;  
g38     = step(G38)';
g38(1) = [ ];
Nss_6 = 76; 
g41     = step(G41)';
Nss_7 = 175;  
g42     = step(G42)';
Nss_8 = 171;  
g44     = step(G44)';
g44(1) = [ ];
Nss_9 = 12;  
g48     = step(G48)';
g48(1) = [ ];
Nss_10 = 125;  
g54     = step(G54)';
g54(1) = [ ];
Nss_11 = 44;  
g56     = step(G56)';
g56(1) = [ ];
Nss_12 = 20;  
g57     = step(G57)';
g57(1) = [ ];
Nss_13 = 47;  
g66    = step(G66)';
g66(1) = [ ];
Nss_14 = 7;  
g71     = step(G71)';
g71(1) = [ ];
Nss_15 = 81;  
g81     = step(G81)';
g81(1) = [ ];
Nss_16 = 54;  
g82     = step(G82)';
g82(1) = [ ];
Nss_17 = 80;  
g83     = step(G83)';
g83(1) = [ ];
Nss_18 =65;  
g88     = step(G88)';
Nss_19 = 31;  

Nss = [Nss_1 Nss_2 Nss_3  Nss_4 Nss_5 Nss_6 Nss_7 ...
            Nss_8  Nss_9 Nss_10 Nss_11 Nss_12 Nss_13  ...
            Nss_14 Nss_15 Nss_16 Nss_17 Nss_18 Nss_19]; 
        
 %% Vetor de coeficientes de Gdi
gd12  =  step(Gd12)'; 
%gd12(1) = [ ];
Nss_d1 =  137;

gd21  =  step(Gd21)'; 
gd21(1) = [ ];
Nss_d2 =  20;

gd22  =  step(Gd22)'; 
gd22(1) = [ ];
Nss_d3 =  68;

gd42  =  step(Gd42)'; 
Nss_d4 =  137;

gd53  =  step(Gd53)'; 
Nss_d5 =  63;

gd63  =  step(Gd63)'; 
gd63(1) = [ ];
Nss_d6 =  7;

gd71  =  step(Gd71)'; 
gd71(1) = [ ];
Nss_d7 = 115 ;

gd72  =  step(Gd72)'; 
gd72(1) = [ ];
Nss_d8 =  23;

gd81  =  step(Gd81)'; 
gd81(1) = [ ];
Nss_d9 =  69;

gd82  =  step(Gd82)'; 
gd82(1) = [ ];
Nss_d10 =  84;
 
Nss_d = [Nss_d1 Nss_d2 Nss_d3  Nss_d4  Nss_d5 ...
                Nss_d6 Nss_d7 Nss_d8 Nss_d9 Nss_d10]; 

%% Define os parametros de sintonia do Controle Preditivo

%ajuste(funcao transferencia, atraso, periodo de amostragem)

%G11
[N1_g11, N2_g11, N_g11, Nu_g11, t5_g11] = ajuste(G11, G11.InputDelay, TS);

%G12
[N1_g12, N2_g12, N_g12, Nu_g12, t5_g12] = ajuste(G12, G12.InputDelay, TS);

%G18
[N1_g18, N2_g18, N_g18, Nu_g18, t5_g18] = ajuste(G18, G18.InputDelay, TS);

%G25
[N1_g25, N2_g25, N_g25, Nu_g25, t5_g25] = ajuste(G25, G25.InputDelay, TS);

%G35
[N1_g35, N2_g35, N_g35, Nu_g35, t5_g35] = ajuste(G35, G35.InputDelay, TS);

%G38
[N1_g38, N2_g38, N_g38, Nu_g38,  t5_g38] = ajuste(G38, G38.InputDelay, TS);

%G41
[N1_g41, N2_g41, N_g41, Nu_g41, t5_g41] = ajuste(G41, G41.InputDelay, TS);

%G42
[N1_g42, N2_g42, N_g42, Nu_g42, t5_g42] = ajuste(G42, G42.InputDelay,TS);

%G44
[N1_g44, N2_g44, N_g44, Nu_g44, t5_g44] = ajuste(G44, G44.InputDelay,TS);

%G48
[N1_g48, N2_g48, N_g48, Nu_g48, t5_g48] = ajuste(G48, G48.InputDelay, TS);

%G54
[N1_g54, N2_g54, N_g54, Nu_g54, t5_g54] = ajuste(G54, G54.InputDelay, TS);

%G56
[N1_g56, N2_g56, N_g56, Nu_g56, t5_g56] = ajuste(G56, G56.InputDelay, TS);

%G57
[N1_g57, N2_g57, N_g57, Nu_g57, t5_g57] = ajuste(G57, G57.InputDelay, TS);

%G66
[N1_g66, N2_g66, N_g66, Nu_g66, t5_g66] = ajuste(G66, G66.InputDelay, TS);

%G71
[N1_g71, N2_g71, N_g71, Nu_g71, t5_g71] = ajuste(G71, G71.InputDelay, TS);

%G81
[N1_g81, N2_g81, N_g81, Nu_g81, t5_g81] = ajuste(G81, G81.InputDelay, TS);

%G82
[N1_g82, N2_g82, N_g82, Nu_g82, t5_g82] = ajuste(G82, G82.InputDelay, TS);

%G83
[N1_g83, N2_g83, N_g83, Nu_g83, t5_g83] = ajuste(G83, G83.InputDelay, TS);

%G88
[N1_g88, N2_g88, N_g88, Nu_g88, t5_g88] = ajuste(G88, G88.InputDelay, TS);


NU = [Nu_g11 Nu_g12 0 0 0 0 0 Nu_g18;
          0 0 0 0 Nu_g25 0 0 0 ;
          0 0 0 0 Nu_g35 0 0 Nu_g38;
          Nu_g41 Nu_g42 0 Nu_g44 0 0 0 Nu_g48;
          0 0 0  Nu_g54 0 Nu_g56 Nu_g57 0;
          0 0 0 0 0 Nu_g66 0 0;
          Nu_g71 0 0 0 0 0 0 0;
          Nu_g81 Nu_g82 Nu_g83 0 0 0 0 Nu_g88];

%% Levantamento das matrizes dinamicas G

[GG11] = matrizDin(N_g11, Nu_g11, G11.InputDelay, g11) ;
[GG41] = matrizDin(N_g41, Nu_g41, G41.InputDelay, g41) ;
[GG71] = matrizDin(N_g71, Nu_g71, G71.InputDelay, g71) ;
[GG81] = matrizDin(N_g81,Nu_g81, G81.InputDelay, g81);
  
[GG12] = matrizDin(N_g12, Nu_g12, G12.InputDelay, g12) ;
[GG42] = matrizDin(N_g42, Nu_g42, G42.InputDelay, g42) ;
[GG82] = matrizDin(N_g82, Nu_g82, G82.InputDelay, g82) ;
 
[GG83] = matrizDin(N_g83, Nu_g83, G83.InputDelay, g83) ;
  
[GG44] = matrizDin(N_g44, Nu_g44, G44.InputDelay, g44);
[GG54] = matrizDin(N_g54, Nu_g54, G54.InputDelay, g54) ;

[GG25] = matrizDin(N_g25, Nu_g25, G25.InputDelay, g25) ;
[GG35] = matrizDin(N_g35, Nu_g35, G35.InputDelay, g35) ;

[GG56] = matrizDin(N_g56, Nu_g56, G56.InputDelay, g56) ;
[GG66] = matrizDin(N_g66, Nu_g66, G66.InputDelay, g66) ;

[GG57] = matrizDin(N_g57, Nu_g57, G57.InputDelay, g57) ;

[GG18] = matrizDin(N_g18, Nu_g18, G18.InputDelay, g18) ;
[GG38] = matrizDin(N_g38, Nu_g38, G38.InputDelay, g38) ;
[GG48] = matrizDin(N_g48, Nu_g48, G48.InputDelay, g48) ;
[GG88] = matrizDin(N_g88, Nu_g88, G88.InputDelay, g88);

%% Redimensionando as matrizes em numero de linhas e colunas

GLen =  [
          [size(GG11,1) size(GG12,1) 0 0 0 0    0    size(GG18,1)]; 
          [0  0    0   0    size(GG25,1)   0   0   0];
          [0  0  0  0    size(GG35,1)   0  0   size(GG38,1)];
          [size(GG41,1) size(GG42,1) 0 size(GG44,1)  0 0 0 size(GG48,1)];
          [0 0 0 size(GG54,1) 0    size(GG56,1)  size(GG57,1)   0];
          [0    0   0   0   0    size(GG66,1)   0   0];
          [size(GG71,1)   0   0   0  0    0    0   0];
          [size(GG81,1) size(GG82,1) size(GG83,1) 0  0 0  0  size(GG88,1)]
          ];
  
GG12 = reshapeMatrix(GLen, GG12, 1);
GG18 = reshapeMatrix(GLen, GG18, 1);
GG35 = reshapeMatrix(GLen, GG35, 3);
GG42 = reshapeMatrix(GLen, GG42, 4);
GG44 = reshapeMatrix(GLen, GG44, 4);
GG48 = reshapeMatrix(GLen, GG48, 4);
GG54 = reshapeMatrix(GLen, GG54, 5);
GG56 = reshapeMatrix(GLen, GG56, 5); 
GG83 = reshapeMatrix(GLen, GG83, 8);
GG88 = reshapeMatrix(GLen, GG88, 8);

GLenDepois =  [
          [size(GG11,1) size(GG12,1) 0 0 0 0    0    size(GG18,1)]; 
          [0  0    0   0    size(GG25,1)   0   0   0];
          [0  0  0  0    size(GG35,1)   0  0   size(GG38,1)];
          [size(GG41,1) size(GG42,1) 0 size(GG44,1)  0 0 0 size(GG48,1)];
          [0 0 0 size(GG54,1) 0    size(GG56,1)  size(GG57,1)   0];
          [0    0   0   0   0    size(GG66,1)   0   0];
          [size(GG71,1)   0   0   0  0    0    0   0];
          [size(GG81,1) size(GG82,1) size(GG83,1) 0  0 0  0  size(GG88,1)]
          ];

GG41 = reshapeCol(NU, GG41, 1);
GG71 = reshapeCol(NU, GG71, 1); 
GG81 = reshapeCol(NU, GG81, 1);
GG42 = reshapeCol(NU, GG42, 2);
GG82 = reshapeCol(NU, GG82, 2); 
GG44 = reshapeCol(NU, GG44, 4);
GG25 = reshapeCol(NU, GG25, 5);
GG18 = reshapeCol(NU, GG18, 8);
GG48 = reshapeCol(NU, GG48, 8);
GG88 = reshapeCol(NU, GG88, 8);

%% Matriz DInamica G

G=   [ GG11                GG12                zeros(126,12)  zeros(126,4)  zeros(126,19)  zeros(126,2)  zeros(126,8)   GG18; 
          zeros(17,45)     zeros(17,26)     zeros(17,12)    zeros(17,4)   GG25                zeros(17,2)     zeros(17,8)     zeros(17,19);
          zeros(53,45)     zeros(53,26)     zeros(53,12)    zeros(53,4)   GG35                zeros(53,2)     zeros(53,8)     GG38;
          GG41                GG42                zeros(62,12)    GG44            zeros(62,19)     zeros(62,2)     zeros(62,8)     GG48;
          zeros(22,45)     zeros(22,26)    zeros(22,12)    GG54             zeros(22,19)     GG56              GG57              zeros(22,19);
          zeros(3,45)       zeros(3,26)      zeros(3,12)       zeros(3,4)      zeros(3,19)       GG66              zeros(3,8)       zeros(3,19)
          GG71                zeros(40,26)    zeros(40,12)    zeros(40,4)    zeros(40,19)     zeros(40,2)     zeros(40,8)     zeros(40,19);
          GG81                GG82               GG83                zeros(40,4)   zeros(40,19)     zeros(40,2)     zeros(40,8)     GG88
       ];

%% Calculo da matriz K (primeira linha)

%horizontes de controle
Nu   = [45 26 12 4 19 2 8 19]; 
%horizontes de predição 
N   = [126 17   53 62  22  3 40 40]';

%normalizando lambdas
%normLambda(ss, N, lambda_escolhido)

% lambda1 = normLambda(ss, du(1), 2);  %peso maior para y1
% lambda2 = normLambda(ss, du(2), 1);
% lambda3 = normLambda(ss, du(3), 1); 
% lambda4 = normLambda(ss, du(4), 1);
% lambda5 = normLambda(ss, du(5), 1); 
% lambda6 = normLambda(ss, du(6), 1); 
% lambda7 = normLambda(ss, du(7), 2); 
% lambda8 = normLambda(ss, du(8), 2); 

% delta1 = normDelta(ss, N(1), 200);  %peso maior para y1
% delta2 = normDelta(ss, N(2), 0.5);
% delta3 = normDelta(ss, N(3), 1); 
% delta4 = normDelta(ss, N(4), 1);
% delta5 = normDelta(ss, N(5), 1); 
% delta6 = normDelta(ss, N(6), 1); 
% delta7 = normDelta(ss, N(7), 200); %peso maior para y7 
% delta8 = normDelta(ss, N(8), 200); %peso maior para y8

lambda1 = 500; lambda2 = 100; lambda3 = 100; lambda4 = 100; lambda5 = 100; lambda6 = 100; lambda7 = 100; lambda8 = 100;

delta1     =  10; delta2 = 100;     delta3 = 100;     delta4 = 1;     delta5 = 10;      delta6=1;       delta7=100;       delta8=100;

Qu = [  lambda1*eye(Nu(1))            zeros(45,26)                     zeros(45,12)                     zeros(45,4)                           zeros(45,19)                         zeros(45,2)                         zeros(45,8)                            zeros(45,19);
            zeros(26,45)                         lambda2*eye(Nu(2))         zeros(26,12)                     zeros(26,4)                           zeros(26,19)                         zeros(26,2)                         zeros(26,8)                            zeros(26,19);
            zeros(12,45)                         zeros(12,26)                      lambda3*eye(Nu(3))        zeros(12,4)                           zeros(12,19)                         zeros(12,2)                         zeros(12,8)                            zeros(12,19);
            zeros(4,45)                           zeros(4,26)                        zeros(4,12)                       lambda4*eye(Nu(4))           zeros(4,19)                            zeros(4,2)                            zeros(4,8)                              zeros(4,19);
            zeros(19,45)                         zeros(19,26)                      zeros(19,12)                    zeros(19,4)                          lambda5*eye(Nu(5))             zeros(19,2)                          zeros(19,8)                            zeros(19,19);
            zeros(2,45)                           zeros(2,26)                        zeros(2,12)                       zeros(2,4)                            zeros(2,19)                           lambda6*eye(Nu(6))            zeros(2,8)                              zeros(2,19)
            zeros(8,45)                           zeros(8,26)                        zeros(8,12)                       zeros(8,4)                            zeros(8,19)                           zeros(8,2)                            lambda7*eye(Nu(7))             zeros(8,19);
            zeros(19,45)                         zeros(19,26)                      zeros(19,12)                    zeros(19,4)                          zeros(19,19)                         zeros(19,2)                          zeros(19,8)                            lambda8*eye(Nu(8));
            ];
        
Qy = [  delta1*eye(N(1))          zeros(126,17)                zeros(126,53)              zeros(126,62)              zeros(126,22)              zeros(126,3)             zeros(126,40)               zeros(126,40);
            zeros(17,126)               delta2*eye(N(2))           zeros(17,53)                zeros(17,62)                 zeros(17,22)                zeros(17,3)               zeros(17,40)                 zeros(17,40);
            zeros(53,126)                zeros(53,17)                  delta3*eye(N(3))        zeros(53,62)                 zeros(53,22)                zeros(53,3)               zeros(53,40)                 zeros(53,40);
            zeros(62,126)                zeros(62,17)                  zeros(62,53)                delta4*eye(N(4))         zeros(62,22)                zeros(62,3)               zeros(62,40)                 zeros(62,40);
            zeros(22,126)                zeros(22,17)                  zeros(22,53)                zeros(22,62)                 delta5*eye(N(5))        zeros(22,3)               zeros(22,40)                 zeros(22,40);
            zeros(3,126)                  zeros(3,17)                     zeros(3,53)                  zeros(3,62)                   zeros(3,22)                delta6*eye(N(6))      zeros(3,40)                   zeros(3,40)
            zeros(40,126)                zeros(40,17)                   zeros(40,53)               zeros(40,62)                 zeros(40,22)              zeros(40,3)              delta7*eye(N(7))           zeros(40,40);
            zeros(40,126)                zeros(40,17)                   zeros(40,53)               zeros(40,62)                 zeros(40,22)              zeros(40,3)              zeros(40,40)                 delta8*eye(N(8));
            ];

K = inv(G'*Qy*G + Qu)*G'*Qy;

%% Calculando resposta livre

%vetor das matrizes Gs
Gs = [G11 G12 G18 G25 G35 G38 G41 G42 G44 G48 ...
            G54 G56 G57 G66 G71 G81 G82 G83 G88];    
        
N= [126 17 53 62 22 3 40 40];

H_g11 = calculoH(Nss,1, g11, N(1));
H_g12 = calculoH(Nss,2, g12, N(1));
H_g18 = calculoH(Nss,3, g18, N(1));
H_g25 = calculoH(Nss,4, g25, N(2));
H_g35 = calculoH(Nss,5, g35, N(3));
H_g38 = calculoH(Nss,6, g38, N(3));
H_g41 = calculoH(Nss,7, g41, N(4));
H_g42 = calculoH(Nss,8, g42, N(4));
H_g44 = calculoH(Nss,9, g44, N(4));
H_g48 = calculoH(Nss,10, g48,N(4));
H_g54 = calculoH(Nss,11, g54, N(5));
H_g56 = calculoH(Nss,12, g56, N(5));
H_g57 = calculoH(Nss,13, g57, N(5));
H_g66 = calculoH(Nss,14, g66, N(6));
H_g71 = calculoH(Nss,15, g71, N(7));
H_g81 = calculoH(Nss,16, g81, N(8));
H_g82 = calculoH(Nss,17, g82, N(8));
H_g83 = calculoH(Nss,18, g83, N(8));
H_g88 = calculoH(Nss,19, g88, N(8));

sizeH = [size(H_g11,2) size(H_g12,2) 0 0 0 0 0 size(H_g18,2);
               0 0 0 0 size(H_g25,2) 0 0 0;
               0 0 0 0 size(H_g35,2) 0 0 size(H_g38,2);
               size(H_g41,2) size(H_g42,2) 0 size(H_g44,2) 0 0 0 size(H_g48,2);
               0 0 0 size(H_g54,2) 0 size(H_g56,2) size(H_g57,2) 0;
               0 0 0 0 0 size(H_g66,2) 0 0;
               size(H_g71,2) 0 0 0 0 0 0 0;
               size(H_g81,2) size(H_g82,2) size(H_g83,2) 0 0 0 0 size(H_g88,2)];

%% adequando o numero de colunas 

H_g41 = [H_g41 zeros(size(H_g41,1),  max(sizeH(:,1)) - size(H_g41,2))];
H_g71 = [H_g71 zeros(size(H_g71,1),  max(sizeH(:,1)) - size(H_g71,2))];
H_g81 = [H_g81 zeros(size(H_g81,1),  max(sizeH(:,1)) - size(H_g81,2))];
H_g42 = [H_g42 zeros(size(H_g42,1),  max(sizeH(:,2)) - size(H_g42,2))];
H_g82 = [H_g82 zeros(size(H_g82,1),  max(sizeH(:,2)) - size(H_g82,2))];
H_g44 = [H_g44 zeros(size(H_g44,1),  max(sizeH(:,4)) - size(H_g44,2))];
H_g25 = [H_g25 zeros(size(H_g25,1),  max(sizeH(:,5)) - size(H_g25,2))];
H_g66 = [H_g66 zeros(size(H_g66,1),  max(sizeH(:,6)) - size(H_g66,2))];
H_g18 = [H_g18 zeros(size(H_g18,1),  max(sizeH(:,8)) - size(H_g18,2))];
H_g38 = [H_g38 zeros(size(H_g38,1),  max(sizeH(:,8)) - size(H_g38,2))];
H_g88 = [H_g88 zeros(size(H_g88,1),  max(sizeH(:,8)) - size(H_g88,2))];
           
H=   [ H_g11              H_g12              zeros(126,65)  zeros(126,44)  zeros(126,85)  zeros(126,20)  zeros(126,47) H_g18; 
        zeros(17,307)   zeros(17,197)  zeros(17,65)    zeros(17,44)    H_g25               zeros(17,20)    zeros(17,47)   zeros(17,125);
        zeros(53,307)   zeros(53,197)  zeros(53,65)    zeros(53,44)    H_g35               zeros(53,20)    zeros(53,47)   H_g38;
        H_g41                H_g42              zeros(62,65)    H_g44              zeros(62,85)     zeros(62,20)    zeros(62,47)   H_g48;
        zeros(22,307)   zeros(22,197)  zeros(22,65)    H_g54              zeros(22,85)     H_g56               H_g57             zeros(22,125);
        zeros(3,307)     zeros(3,197)    zeros(3,65)      zeros(3,44)       zeros(3,85)       H_g66               zeros(3,47)     zeros(3,125)
        H_g71               zeros(40,197)  zeros(40,65)    zeros(40,44)    zeros(40,85)     zeros(40,20)     zeros(40,47)    zeros(40,125);
        H_g81               H_g82               H_g83              zeros(40,44)    zeros(40,85)     zeros(40,20)     zeros(40,47)   H_g88
        ];
  
 %% Calculando H para as perturbacoes 
 Hd_g12 = calculoH(Nss_d,1, gd12, N(1));
 Hd_g21 = calculoH(Nss_d,2, gd21, N(2));
 Hd_g22 = calculoH(Nss_d,3, gd22, N(2));
 Hd_g42 = calculoH(Nss_d,4, gd42, N(4));
 Hd_g53 = calculoH(Nss_d,5, gd53, N(5));
 Hd_g63 = calculoH(Nss_d,6, gd63, N(6));
 Hd_g71 = calculoH(Nss_d,7, gd71, N(7));
 Hd_g72 = calculoH(Nss_d,8, gd72, N(7));
 Hd_g81 = calculoH(Nss_d,9, gd81,   N(8));
 Hd_g82 = calculoH(Nss_d,10, gd82, N(8));
 
 sizeHd = [
           0        size(Hd_g12,2)    0;
          size(Hd_g21,2)  size(Hd_g22,2)    0;
           0         0            0;
           0         size(Hd_g42,2)     0;
           0         0            size(Hd_g53,2);
           0         0            size(Hd_g63,2);
           size(Hd_g71,2) size(Hd_g72,2)     0;
           size(Hd_g81,2) size(Hd_g82,2)     0
           ];
 
 %ajustando colunas
Hd_g21 = [Hd_g21 zeros(size(Hd_g21,1),  max(sizeHd(:,1)) - size(Hd_g21,2))];
Hd_g81 = [Hd_g81 zeros(size(Hd_g81,1),  max(sizeHd(:,1)) - size(Hd_g81,2))];
Hd_g22 = [Hd_g22 zeros(size(Hd_g22,1),  max(sizeHd(:,2)) - size(Hd_g22,2))];
Hd_g72 = [Hd_g72 zeros(size(Hd_g72,1),  max(sizeHd(:,2)) - size(Hd_g72,2))];
Hd_g82 = [Hd_g82 zeros(size(Hd_g82,1),  max(sizeHd(:,2)) - size(Hd_g82,2))];
Hd_g63 = [Hd_g63 zeros(size(Hd_g63,1),  max(sizeHd(:,3)) - size(Hd_g63,2))]; 

Nssd = [115 137 63];

Hq = [
           zeros(N(1),Nssd(1))    Hd_g12                     zeros(N(1),Nssd(3));
           Hd_g21                        Hd_g22                     zeros(N(2),Nssd(3));
           zeros(N(3),Nssd(1))    zeros(N(3),Nssd(2))  zeros(N(3),Nssd(3));
           zeros(N(4),Nssd(1))    Hd_g42                      zeros(N(4),Nssd(3));
           zeros(N(5),Nssd(1))    zeros(N(5),Nssd(2))  Hd_g53;
           zeros(N(6),Nssd(1))    zeros(N(6),Nssd(2))  Hd_g63;
           Hd_g71                        Hd_g72                     zeros(N(7),Nssd(3));
           Hd_g81                        Hd_g82                     zeros(N(8),Nssd(3))
           ];
       
 %%%%%%%%%SIMULACAO%%%%%%%%%%%%%%%%
%% Condicoes iniciais

%iteracoes
iters = 5000;

%matriz diagonal I
I = [ones(126,1) zeros(126,1) zeros(126,1) zeros(126,1) zeros(126,1) zeros(126,1) zeros(126,1) zeros(126,1);
      zeros(17,1)  ones(17,1) zeros(17,1) zeros(17,1) zeros(17,1) zeros(17,1) zeros(17,1) zeros(17,1);
      zeros(53,1)   zeros(53,1) ones(53,1) zeros(53,1) zeros(53,1) zeros(53,1) zeros(53,1) zeros(53,1);
      zeros(62,1)   zeros(62,1) zeros(62,1) ones(62,1) zeros(62,1) zeros(62,1) zeros(62,1) zeros(62,1);
      zeros(22,1)  zeros(22,1)  zeros(22,1) zeros(22,1) ones(22,1) zeros(22,1) zeros(22,1) zeros(22,1);
      zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)  ones(3,1) zeros(3,1)  zeros(3,1) ;
     zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) ones(40,1) zeros(40,1);
      zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) zeros(40,1) ones(40,1)
      ]; 

%horizontes de modelo e  inicial
Nss = [307 197 65 44 85 20 47 125];
N1   = [1 1 1 1 1 1 1 1]';

%vetores delta u passados da resposta livre
%dimensao 1xNss
du_pass1 = zeros(Nss(1),1);
du_pass2 = zeros(Nss(2),1);
du_pass3 = zeros(Nss(3),1);
du_pass4 = zeros(Nss(4),1);
du_pass5 = zeros(Nss(5),1);
du_pass6 = zeros(Nss(6),1);
du_pass7 = zeros(Nss(7),1);
du_pass8 = zeros(Nss(8),1);
du_passados = [du_pass1; du_pass2; du_pass3; du_pass4; du_pass5; du_pass6; du_pass7; du_pass8];


%vetores delta q passados da resposta livre
dq_pass1 = zeros(Nssd(1),1);
dq_pass2 = zeros(Nssd(2),1);
dq_pass3 = zeros(Nssd(3),1);
dq_passados = [dq_pass1; dq_pass2; dq_pass3];

%vetor do sinal de controle
u = zeros(8, iters);

%inicializando variavel de incremento de controle
%dimensao 1xNu,l + 1xNu,2 + ... 1xNu,8
deltaU = zeros(sum(Nu),1);

%vetores para guardar saidas 
out = [ ]; 
y = zeros(sum(N), 1);
ysim = zeros(8,iters);

%vetor do tempo
t = 1:iters;

%vetores das referencias
w1 =  [zeros(1,100)  8*ones(1,4900)];
w2 =  [zeros(1,200) 7.5*ones(1,4800)];
w3 =  [zeros(1,300) 7*ones(1, 4700)];
w4 =  [zeros(1,400) 6.5*ones(1,4600)];
w5 =  [zeros(1,500) 6*ones(1,4500)];
w6 =  [zeros(1,600) 5.5*ones(1,4400)];
w7 =  [zeros(1,700) 5*ones(1,4300)];
w8 =  [zeros(1,800) 4.5*ones(1,4200)];

w = [w1;w2;w3;w4;w5;w6;w7;w8];
 
%% DMC Loop

for k=2:iters  

inputs = u';
out = lsim(Gu, inputs, t);
ysim(:,k) = out(k,:)';

y = G*deltaU + H*du_passados + I*ysim(:,k) + Hq*dq_passados;

 %atualizacao deltaUs passados
 du_passados = [deltaU(1,1);
                            du_passados(1:Nss(1)-1,1);
                            deltaU(Nu(1)+1,1);
                            du_passados(Nss(1)+1:Nss(1)+Nss(2)-1,1);
                            deltaU(Nu(1)+Nu(2)+1,1);
                            du_passados(Nss(1)+Nss(2)+1:Nss(1)+Nss(2)+Nss(3)-1,1);
                            deltaU(Nu(1)+Nu(2)+Nu(3)+1,1);
                            du_passados(Nss(1)+Nss(2)+Nss(3)+1:Nss(1)+Nss(2)+Nss(3)+Nss(4)-1,1);
                            deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+1,1);
                            du_passados(Nss(1)+Nss(2)+Nss(3)+Nss(4)+1:Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)-1,1);
                            deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+1,1);
                            du_passados(Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+1:Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+Nss(6)-1,1);
                            deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+Nu(6)+1,1);
                            du_passados(Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+Nss(6)+1:Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+Nss(6)+Nss(7)-1,1);
                            deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+Nu(6)+Nu(7)+1,1);
                            du_passados(Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+Nss(6)+Nss(7)+1:Nss(1)+Nss(2)+Nss(3)+Nss(4)+Nss(5)+Nss(6)+Nss(7)+Nss(8)-1,1);
                            ];
                        
%atualizacao deltaQs passados
%  dq_passados = [dq_passados(1:Nssd(1),1);
%                              dq_passados(Nssd(1)+1:Nssd(1)+Nssd(2),1);
%                              dq_passados(Nssd(1)+Nssd(2)+1:Nssd(1)+Nssd(2)+Nssd(3),1);
%                              ];
 
%set point  
%considera ref futura cte
ref1 =w1(k)*ones(1,N(1))';
ref2 =w2(k)*ones(1,N(2))';
ref3 =w3(k)*ones(1,N(3))';
ref4 =w4(k)*ones(1,N(4))';
ref5 =w5(k)*ones(1,N(5))';
ref6 =w6(k)*ones(1,N(6))';
ref7 =w7(k)*ones(1,N(7))';
ref8 =w8(k)*ones(1,N(8))';

ref = [ref1 ; ref2 ;  ref3 ; ref4 ; ref5 ; ref6; ref7; ref8];
        
%calcula o controle 
deltaU = K*(ref-y);

%atualiza controle
u(1,k) = u(1, k-1) + deltaU(1,1);
u(2,k) = u(2, k-1) + deltaU(Nu(1)+1,1);
u(3,k) = u(3,k-1) + deltaU(Nu(1)+Nu(2)+1,1);
u(4,k) = u(4,k-1) + deltaU(Nu(1)+Nu(2)+Nu(3)+1,1);
u(5,k) = u(5,k-1) + deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+1,1);
u(6,k) = u(6,k-1) + deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+1,1);
u(7,k) = u(7,k-1) + deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+Nu(6)+1,1);
u(8,k) = u(8,k-1) + deltaU(Nu(1)+Nu(2)+Nu(3)+Nu(4)+Nu(5)+Nu(6)+Nu(7)+1,1);

end


%%  Graficos

%saidas e referencias

figure(1)

%saida 1, ref 1
subplot(8,1,1)
plot(t,out(:,1), 'r')
hold on
plot(t,w(1,:), 'r')
legend('saida1','ref1')

%saida 2, ref 2
subplot(8,1,2)
plot(t,out(:,2), 'g')
hold on
plot(t,w(2,:), 'g')
legend('saida2','ref2')

%saida 3, ref 3
subplot(8,1,3)
plot(t,out(:,3), 'b')
hold on
plot(t,w(3,:), 'b')
legend('saida3','ref3')

%saida 4, ref 4
subplot(8,1,4)
plot(t,out(:,4), 'k')
hold on
plot(t,w(4,:),'k')
legend('saida4','ref4')

%saida 5, referencia 5
subplot(8,1,5)
plot(t,out(:,5), 'm')
hold on
plot(t,w(5,:),'m')
legend('saida5','ref5')

%saida 6 referencia 6
subplot(8,1,6)
plot(t,out(:,6), 'y')
hold on
plot(t,w(6,:),'y')
legend('saida6','ref6')

%saida 7, referencia 7
subplot(8,1,7)
plot(t,out(:,7), 'c')
hold on
plot(t,w(7,:),'c')
legend('saida7','ref7')

%saida 8, referencia 8
subplot(8,1,8)
plot(t,out(:,8), 'k')
hold on
plot(t,w(8,:),'k')
legend('saida8','ref8')

%sinal de controle
figure(2)

%u1
subplot(8,1,1);
plot(t,u(1,:),'r')
legend('u1')
%u2
subplot(8,1,2);
plot(t,u(2,:),'g')
legend('u2')
%u3
subplot(8,1,3);
plot(t,u(3,:),'b')
legend('u3')
%u4
subplot(8,1,4);
plot(t,u(4,:),'k')
legend('u4')
%u5
subplot(8,1,5);
plot(t,u(5,:),'m')
legend('u5')
%u6
subplot(8,1,6);
plot(t,u(6,:),'y')
legend('u6')
%u7
subplot(8,1,7);
plot(t,u(7,:),'c')
legend('u7')
%u8
subplot(8,1,8);
plot(t,u(8,:),'k')
legend('u8')


% subplot(3,1,3)
% p = [zeros(1,399) 0.5*ones(1,iters)];
% plot(p(1:500), 'g');
% legend('perturbacao')

