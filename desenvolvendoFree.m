clc 
clear all

H = [ ] ;

format shortg

N2 =    5;
Nss = 10;

%h = [ 0.1 0.4 0.6 0.7 0.85 0.95 1 1.2 1.5 1.8 2 2.75 3 4.5 5 6.6];
%h = [ 0.1 0.4 0.6 0.7  0.85 0.95];
h = [ 0.1 0.4 0.6];

z = 1 + N2 + Nss;

if(z > length(h))
    sprintf("Valor %d", z)
    h = [h, h(end)*ones(1,z-length(h))];
else
    N2 = 5;
    Nss = 10;
end

for i = 1:Nss
    g(i) = h(1+i) - h(i);
end
 H(1,:) = g;

for k = 1:N2-1
    for j = 1:Nss
        g(j) = h(1+j+k) - h(j);
        H(k+1,:) = g; 
    end
end
