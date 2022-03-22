function [ K1 ] = calculandoK( delta, G, lambda, N, Nu)

K =inv(G'*delta*eye(N)*G + lambda*eye(Nu))*G';
K1=  K(1,:);
    
end

