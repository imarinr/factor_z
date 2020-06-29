
% Correlación para obtener el factor de compresibilidad de los gases reales
% por Dranchuk y Abou-Kassem
% Tpc es la temperatura pseudocrítica en °R,
% T es la temperatura a condiciones de yacimiento en °R,
% Tpr es la temperatura pseudoreducida en °R,
% Ppc es la presión pseudocrítica en lb/pg2 absolutas,
% P es la presión en lb/pg2 absolutas,
% Ppr es la presión pseudoreducida en lb/pg2 absolutas,
% z es el factor de compresibilidad su valor inicial es igual a 1
% Rhor  es la densidad pseudoreducida en lbm/pie3.

% MODO DE USO:
% En la ventana de MATLAB ubiquese en la carpeta que contiene este archivo
% Escriba:
% >> factorZ(temperatura_pseudoreducida, presion_pseudoreducida)
% y presione enter
function z=factorZ(Tpr, Ppr)
    % Definicion de constantes
    a1=0.3265;
    a2=-1.07;
    a3=-0.5339;
    a4=0.01569;
    a5=-0.05165;
    a6=0.5475;
    a7=-0.7361;
    a8=0.1844;
    a9=0.1056;
    a10=0.6134;
    a11=0.721;

    % Temperatura y presion pseudocriticas
    %Tpc = 168 + 325*densRelativa -  12.5*densRelativa^2;
    %Ppc = 677 + 15*densRelativa - 37.5*densRelativa^2;

    % Temperatura y presion pseudoreducidas
    %Tpr=temp/Tpc;
    %Ppr=presion/Ppc;

    %Si 1 < Tpr <= 3 y  0.2 < Ppr < 30 continua, de lo contrario no aplica
    if and(Tpr > 1, Tpr <= 3)
        if and(Ppr >= 0.2, Ppr < 30)
            % El valor inical de z es la unidad
            z=1;

            Rhor=0.27*Ppr/(z*Tpr);

             c1=a1+(a2/Tpr)+(a3/(Tpr^3))+(a4/(Tpr^4))+(a5/(Tpr^5));
             c2=a6+(a7/Tpr)+(a8/(Tpr^2));
             c3=a9*((a7/Tpr)+(a8/(Tpr^2)));

             c4=a10*(1+a11*(Rhor^2))*((Rhor^2)/(Tpr^3))*exp(-a11*(Rhor^2));

             % aqui empieza el ciclo ya que comenza en 1 y de ahi va el ciclio
             % evaluando en rho que es la densidad relativa y en c4 que es
             % funcion de la presion y de rhor tambien se evalua la derivada
             % tolerancia indica la precision del resultado maxIter limita a
             % cien iteraciones el metodo, esto para evitar un bucle infinito
             tolerancia = 0.0001;
             maxIter = 100;
             for n = 1:maxIter
                 Fz = z - (1 + (c1 * Rhor) + (c2 * Rhor^2) - (c3 * Rhor^5) + c4);
                 dFz = 1 + (c1 * (Rhor / z)) + ((2 * c2) * (Rhor^2 / z)) - ((5 * c3) * (Rhor^5 / z)) + ((2*a10* Rhor^2 / (Tpr^3 * z)) * (1 + (a11 * Rhor^2) - (a11 * Rhor^2)^2) * exp(-a11*Rhor^2));
                 znew = z -(Fz / dFz);
                 Rhor=0.27*Ppr/(znew*Tpr);
                 c4=a10*(1+a11*(Rhor^2))*((Rhor^2)/(Tpr^3))*exp(-a11*(Rhor^2));
                 dif = abs(z - znew);

                 if or(dif < tolerancia, Fz == 0)
                     break;
                 end

                 z = znew;
             end
        end
    else
        fprintf('El método no aplica\n')
    end
