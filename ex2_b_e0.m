%Definições de parâmetros iniciais
m = 500;
k_y = 20000;
c_y = 700;

T = 15; %Tempo de simulação
h = 0.0001;
V = 3; 
  
t = 0:h:T; %Vetor de tempo

x1 = zeros(1, length(t)); % y
x2 = zeros(1, length(t)); % dy/dt
x2_dot = zeros(1, length(t)); % d2y/dt2

% Condições iniciais
x1(1)=0;
x2(1)=0;


%Resolução da EDO
for i = 1:length(t)-1
    
    [z, z_dot] = lombada(t(i), V);
    
    %Obtenção dos k_i's intermediários
    k_1 = zeros(1, 4);
    k_2 = zeros(1, 4);

    k_1(1) = x2(i);
    k_2(1) = (-1/m)*(k_y*(x1(i)-z)+c_y*(x2(i)-z_dot));

    for v = 2:4
        x1_aux = x1(i) + h*k_1(v-1);
        x2_aux = x2(i) + h*k_2(v-1);
        x2_dot_aux = (-1/m)*(k_y*(x1_aux-z)+c_y*(x2_aux-z_dot));
        k_1(v) = x2_aux;
        k_2(v) = x2_dot_aux;
    end

f1 = (k_1(1)+2*k_1(2)+2*k_1(3)+k_1(4))/6;
f2 = (k_2(1)+2*k_2(2)+2*k_2(3)+k_2(4))/6;

x1(i+1)=x1(i) + h*f1;
x2(i+1)=x2(i) + h*f2;
x2_dot(i) = f2;

end

% Plot
figure('Position', [100, 100, 1000, 400])
plot(t, x1, 'LineWidth', 2, 'DisplayName', 'Sem absorvedor')
hold on; title("Deslocamento"); ylabel('y (m)'); legend
plot(t, Y(1,:), 'LineWidth', 2, 'DisplayName', 'Com absorvedor ')

function [z,z_dot]=lombada(t, V)
    
    T = 2/V;
    
    if t < floor(2/T)*T
        z = 0;
        z_dot = 0;
    elseif t < (floor(2/T)+10)*T
        z = 0.125*(1-cos(2*pi*t/T));
        z_dot = (0.25*pi/T)*sin(2*pi*t/T);
    else
        z = 0;
        z_dot = 0;
    end

end


