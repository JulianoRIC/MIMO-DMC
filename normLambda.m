function [ lambda_norm ] = normLambda( deltaUs, N, lambda)


 %normalizando lambda
 lambda_norm = lambda/(  (max(deltaUs)^2)*(N));

end

