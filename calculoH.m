function [ H ] = calculoH( nss, nro, h, n )

H = [ ] ;

N =    n;
Nss = nss(nro);
z = 1 + N + Nss;

if(z > length(h))
    sprintf("Valor %d", z);
    h = [h, h(end)*ones(1,z-length(h))];
% else
%     N = n;
%     Nss = nss;
end

for i = 1:Nss
    g(i) = h(1+i) - h(i);
end
 H(1,:) = g;

for k = 1:N-1
    for j = 1:Nss
        g(j) = h(1+j+k) - h(j);
        H(k+1,:) = g; 
    end
end
