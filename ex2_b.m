clc, clear

T = 15; %Tempo de simulação
h = 0.0001; %Passo

t = 0:h:T; %Vetor de tempo

x1 = zeros(1, length(t)); % y
x2 = zeros(1, length(t)); % dy/dt
x2_dot = zeros(1, length(t)); % d2y/dt2
x3 = zeros(1, length(t)); % u
x4 = zeros(1, length(t)); % du/dt
x4_dot = zeros(1, length(t)); % d2u/dt2

Y = [x1; x2; x3; x4];

% Condições iniciais
x1(1)=0;
x2(1)=0;
x3(1)=0;
x4(1)=0;


epsilon_vals = 0.01:0.045:0.19;

figure('Position', [100, 100, 700, 600])

for epsilon = epsilon_vals
   
    for i = 1:length(t)-1
    
        K1 = f(t(i), Y(:,i), epsilon);
        K2 = f(t(i)+h/2 , Y(:,i)+(h/2)*K1, epsilon);
        K3 = f(t(i)+h/2 , Y(:,i)+(h/2)*K2, epsilon);
        K4 = f(t(i)+h , Y(:,i)+h*K3, epsilon);

        Y(:, i+1) = Y(:,i) + (h/6)*(K1 + 2*K2 + 2*K3 + K4);
    
    end
    
    subplot(2,1,1)
    plot(t, Y(1,:), 'LineWidth', 2, 'DisplayName', sprintf('\\epsilon = %.3f', epsilon))
    hold on; title("Deslocamento de y"); ylabel('y (m)'); legend

    subplot(2,1,2)
    plot(t, Y(3,:), 'LineWidth', 2, 'DisplayName', sprintf('\\epsilon = %.3f', epsilon))
    hold on; title("Deslocamento de u"); ylabel('u (m)'); xlabel('t (s)'); legend
    
end


function funcao=f(t,Y, epsilon)

    M = 500;
    k_y = 20000;
    c_y = 700;
    k = 10000;
    c = 350;
    L = 0.1;
    m = epsilon * M;
    V = 3; %Velocidade de avanço do veiculo
    
    [z, z_dot] = lombada(t,V);

    funcao = [Y(2);
        (1/M)*((k*Y(3)^3)/(2*L^2) + c*Y(4) - k_y*(Y(1)-z) - c_y*(Y(2) - z_dot));
        Y(4);
        (-1/M)*(((M+m)/m)*((k*Y(3)^3)/(2*L^2) + c*Y(4)) - k_y*(Y(1)-z) - c_y*(Y(2) - z_dot))
        ];

end

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









