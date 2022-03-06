function [ G ] = matrizDin( N, Nu, d, gi)

% Calculo da Matriz G
G = zeros(N,Nu);
G(:,1) = gi(d+1:N+d);         
for i =2:Nu
    for j = 2:N
    G(j,i) = G(j-1,i-1);
    end
end

