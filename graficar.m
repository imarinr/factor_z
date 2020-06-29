% esta funcion se encargara de graficar los valores de z en 1 < Tpr <= 3 y
% 0.2 < Ppr < 30 
tr = 1.1:0.1:3
pr = 0.2:0.1:29.9
z = zeros(length(pr))
iMax = length(tr);
jMax = length(pr);
for i=1:iMax
    for j = 1: jMax
        z(j) = factorZ(tr(i), pr(j));
    end
    plot(pr, z)
    hold on;
end

