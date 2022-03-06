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
G1_2 =  tf([0.01366 -0.0128 -0.01706], [1 -1.915 0.9191], TS, 'InputDelay', 20);
G21   =  tf([0.00236 0.002173], [1 -0.9518 0.1157], TS);
G22   =  tf([4.68 -3.722], [1 -0.8786 3.875e-5], TS, 'InputDelay', 20);
G4_2 =  tf([-0.00611 0.005808 0.007654], [1 -1.914 0.9181], TS, 'InputDelay', 20);
G53   =  tf([0.03646 -0.006099 -0.02366], [1 -1.444 0.4712], TS, 'InputDelay', 6);
G63   =  tf([0.5015 0.2255], [1 -0.3526 0.0961], TS);
G7_1 =  tf([0.0004086 -0.0003754], [1 -1.306 0.3331], TS);
G72   =  tf([0.07873 -0.07549], [1 -1.176 0.2297], TS);
G8_1 =  tf([3.153e-5 -2.93810e-5], [1 -1.9  0.903], TS);
G8_2 =  tf([0 0.008119], [1 -0.956], TS);

% matriz perturbacao
Gw = [
           0      G1_2  0;
           G21 G22    0;
           0       0        0;
           0      G4_2  0;
           0      0     G53;
           0      0     G63;
           G7_1 G72   0;
           G8_1 G8_2 0
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
Nss_6 = 76; %atencao
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

%% Define os parametros de sintonia do Controle Preditivo

%ajuste(funcao transferencia, atraso, lambda, delta, periodo de amostragem)

%G11
[N1_g11, N2_g11, N_g11, Nu_g11, lambda_g11, delta_g11, t5_g11] = ajuste(G11, G11.InputDelay, 1, 1, TS);

%G12
[N1_g12, N2_g12, N_g12, Nu_g12, lambda_g12, delta_g12, t5_g12] = ajuste(G12, G12.InputDelay, 1, 1, TS);

%G18
[N1_g18, N2_g18, N_g18, Nu_g18, lambda_g18, delta_g18, t5_g18] = ajuste(G18, G18.InputDelay, 1, 1, TS);

%G25
[N1_g25, N2_g25, N_g25, Nu_g25, lambda_g25, delta_g25, t5_g25] = ajuste(G25, G25.InputDelay, 1, 1, TS);

%G35
[N1_g35, N2_g35, N_g35, Nu_g35, lambda_g35, delta_g35, t5_g35] = ajuste(G35, G35.InputDelay, 1, 1, TS);

%G38
[N1_g38, N2_g38, N_g38, Nu_g38, lambda_g38, delta_g38, t5_g38] = ajuste(G38, G38.InputDelay, 1, 1, TS);

%G41
[N1_g41, N2_g41, N_g41, Nu_g41, lambda_g41, delta_g41, t5_g41] = ajuste(G41, G41.InputDelay, 1, 1, TS);

%G42
[N1_g42, N2_g42, N_g42, Nu_g42, lambda_g42, delta_g42, t5_g42] = ajuste(G42, G42.InputDelay, 1, 1, TS);

%G44
[N1_g44, N2_g44, N_g44, Nu_g44, lambda_g44, delta_g44, t5_g44] = ajuste(G44, G44.InputDelay, 1, 1, TS);

%G48
[N1_g48, N2_g48, N_g48, Nu_g48, lambda_g48, delta_g48, t5_g48] = ajuste(G48, G48.InputDelay, 1, 1, TS);

%G54
[N1_g54, N2_g54, N_g54, Nu_g54, lambda_g54, delta_g54, t5_g54] = ajuste(G54, G54.InputDelay, 1, 1, TS);

%G56
[N1_g56, N2_g56, N_g56, Nu_g56, lambda_g56, delta_g56, t5_g56] = ajuste(G56, G56.InputDelay, 1, 1, TS);

%G57
[N1_g57, N2_g57, N_g57, Nu_g57, lambda_g57, delta_g57, t5_g57] = ajuste(G57, G57.InputDelay, 1, 1, TS);

%G66
[N1_g66, N2_g66, N_g66, Nu_g66, lambda_g66, delta_g66, t5_g66] = ajuste(G66, G66.InputDelay, 1, 1, TS);

%G71
[N1_g71, N2_g71, N_g71, Nu_g71, lambda_g71, delta_g71, t5_g71] = ajuste(G71, G71.InputDelay, 1, 1, TS);

%G81
[N1_g81, N2_g81, N_g81, Nu_g81, lambda_g81, delta_g81, t5_g81] = ajuste(G81, G81.InputDelay, 1, 1, TS);

%G82
[N1_g82, N2_g82, N_g82, Nu_g82, lambda_g82, delta_g82, t5_g82] = ajuste(G82, G82.InputDelay, 1, 1, TS);

%G83
[N1_g83, N2_g83, N_g83, Nu_g83, lambda_g83, delta_g83, t5_g83] = ajuste(G83, G83.InputDelay, 1, 1, TS);

%G88
[N1_g88, N2_g88, N_g88, Nu_g88, lambda_g88, delta_g88, t5_g88] = ajuste(G88, G88.InputDelay, 1, 1, TS);


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

K1_g11 = calculandoK(delta_g11, GG11, lambda_g11, 126, 45) ;
K1_g41 = calculandoK(delta_g41, GG41, lambda_g41, 62, 45) ;
K1_g71 = calculandoK(delta_g71, GG71, lambda_g71, 40, 45) ;
K1_g81 = calculandoK(delta_g81, GG81, lambda_g81, 40, 45) ;
K1_g12 = calculandoK(delta_g12, GG12, lambda_g12, 126, 26) ; 
K1_g42 = calculandoK(delta_g42, GG42, lambda_g42, 62,   26) ;
K1_g82 = calculandoK(delta_g82, GG82, lambda_g82, 40,   26) ;
K1_g83 = calculandoK(delta_g83, GG83, lambda_g83, 40, 12) ;
K1_g44 = calculandoK(delta_g44, GG44, lambda_g44, 62, 4) ;
K1_g54 = calculandoK(delta_g54, GG54, lambda_g54, 22, 4) ;
K1_g25 = calculandoK(delta_g25, GG25, lambda_g25, 17, 19) ;
K1_g35 = calculandoK(delta_g35, GG35, lambda_g35, 53, 19) ;
K1_g56 = calculandoK(delta_g56, GG56, lambda_g56, 22, 2) ;
K1_g66 = calculandoK(delta_g66, GG66, lambda_g66, 3,   2) ;
K1_g57 = calculandoK(delta_g57, GG57, lambda_g57, 22, 8) ;
K1_g18 = calculandoK(delta_g18, GG18, lambda_g18, 126, 19) ;
K1_g38 = calculandoK(delta_g38, GG38, lambda_g38, 53, 19) ;
K1_g48 = calculandoK(delta_g48, GG48, lambda_g48, 62, 19) ;
K1_g88 = calculandoK(delta_g88, GG88, lambda_g88, 40, 19) ;

% 
% %%adequando colunas
K1_g41 = [K1_g41 zeros(1, 126 - 62)];
K1_g71 = [K1_g71 zeros(1, 126 - 40)];
K1_g81 = [K1_g81 zeros(1, 126 - 40)];
K1_g42 = [K1_g42 zeros(1, 126 - 62)];
K1_g82 = [K1_g82 zeros(1, 126 - 40)];
K1_g54 = [K1_g54 zeros(1, 62- 22)];
K1_g25 = [K1_g25 zeros(1, 53- 17)];
K1_g66 = [K1_g66 zeros(1, 22- 3)];
K1_g38 = [K1_g38 zeros(1, 126- 53)];
K1_g48 = [K1_g48 zeros(1, 126- 62)];
K1_g88 = [K1_g88 zeros(1, 126- 40)];

K1 = [K1_g11             K1_g12            zeros(1,40)      zeros(1,62)      zeros(1,53)     zeros(1,22)   zeros(1,22)    K1_g18;
          zeros(1,126)     zeros(1,126)    zeros(1,40)      zeros(1,62)      K1_g25          zeros(1,22)   zeros(1,22)     zeros(1,126);
          zeros(1,126)     zeros(1,126)    zeros(1,40)      zeros(1,62)      K1_g35          zeros(1,22)   zeros(1,22)     K1_g38;
          K1_g41             K1_g42            zeros(1,40)      K1_g44           zeros(1,53)     zeros(1,22)   zeros(1,22)     K1_g48;
          zeros(1,126)     zeros(1,126)    zeros(1,40)      K1_g54           zeros(1,53)     K1_g56         K1_g57           zeros(1,126);
          zeros(1,126)     zeros(1,126)    zeros(1,40)      zeros(1,62)      zeros(1,53)    K1_g66        zeros(1,22)      zeros(1,126);
          K1_g71             zeros(1,126)    zeros(1,40)      zeros(1,62)      zeros(1,53)    zeros(1,22)   zeros(1,22)     zeros(1,126);
          K1_g81             K1_g82            K1_g83            zeros(1,62)      zeros(1,53)    zeros(1,22)   zeros(1,22)     K1_g88
      ];

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

sizeLHD = [size(H_g11,1) size(H_g12,1) 0 0 0 0 0 size(H_g18,1);
               0 0 0 0 size(H_g25,1) 0 0 0;
               0 0 0 0 size(H_g35,1) 0 0 size(H_g38,1);
               size(H_g41,1) size(H_g42,1) 0 size(H_g44,1) 0 0 0 size(H_g48,1);
               0 0 0 size(H_g54,1) 0 size(H_g56,1) size(H_g57,1) 0;
               0 0 0 0 0 size(H_g66,1) 0 0;
               size(H_g71,1) 0 0 0 0 0 0 0;
               size(H_g81,1) size(H_g82,1) size(H_g83,1) 0 0 0 0 size(H_g88,1)];
           
sizeCHD = [size(H_g11,2) size(H_g12,2) 0 0 0 0 0 size(H_g18,2);
               0 0 0 0 size(H_g25,2) 0 0 0;
               0 0 0 0 size(H_g35,2) 0 0 size(H_g38,2);
               size(H_g41,2) size(H_g42,2) 0 size(H_g44,2) 0 0 0 size(H_g48,2);
               0 0 0 size(H_g54,2) 0 size(H_g56,2) size(H_g57,2) 0;
               0 0 0 0 0 size(H_g66,2) 0 0;
               size(H_g71,2) 0 0 0 0 0 0 0;
               size(H_g81,2) size(H_g82,2) size(H_g83,2) 0 0 0 0 size(H_g88,2)];
           
           
H=   [ H_g11              H_g12              zeros(126,65)  zeros(126,44)  zeros(126,85)  zeros(126,20)  zeros(126,47) H_g18; 
        zeros(17,307)   zeros(17,197)  zeros(17,65)    zeros(17,44)    H_g25               zeros(17,20)    zeros(17,47)   zeros(17,125);
        zeros(53,307)   zeros(53,197)  zeros(53,65)    zeros(53,44)    H_g35               zeros(53,20)    zeros(53,47)   H_g38;
        H_g41                H_g42              zeros(62,65)    H_g44              zeros(62,85)     zeros(62,20)    zeros(62,47)   H_g48;
        zeros(22,307)   zeros(22,197)  zeros(22,65)    H_g54              zeros(22,85)     H_g56               H_g57             zeros(22,125);
        zeros(3,307)     zeros(3,197)    zeros(3,65)      zeros(3,44)       zeros(3,85)       H_g66               zeros(3,47)     zeros(3,125)
        H_g71               zeros(40,197)  zeros(40,65)    zeros(40,44)    zeros(40,85)     zeros(40,20)     zeros(40,47)    zeros(40,125);
        H_g81               H_g82               H_g83              zeros(40,44)    zeros(40,85)     zeros(40,20)     zeros(40,47)   H_g88
        ];

%iteracoes
iters = 500;