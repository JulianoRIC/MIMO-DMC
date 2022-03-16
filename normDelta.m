function [ delta_norm ] = normDelta( ss, N, delta)

 %normalizando lambda
 delta_norm = delta/(  (max(ss)^2)*(N));

end

