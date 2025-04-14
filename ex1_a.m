%Definições de parâmetros iniciais
m = 500;
k_y = 20000;
c_y = 700;

T = 5; %Tempo de simulação
h_vals = [0.01, 0.001, 0.0001];


% A EDO será resolvida sendo passada para o espaço de estados.
% A forma matricial não será usada para facilitar plots e por se tratar de
% poucas variáveis.

for h = h_vals
  
    t = 0:h:T; %Vetor de tempo

    x1 = zeros(1, length(t)); % y
    x2 = zeros(1, length(t)); % dy/dt
    x2_dot = zeros(1, length(t)); % d2y/dt2

    % Condições iniciais
    x1(1)=0;
    x2(1)=0;
    z_dot = 0; % A derivada do degrau é 0 em todo ponto, menos em t=2, que é infinito.

    %Resolução da EDO
    for i = 1:length(t)-1

        f1 = x2(i);
        f2 = (-1/m)*(k_y*(x1(i)-z(t(i)))+c_y*(x2(i)-z_dot));

        x1(i+1)=x1(i) + h*f1;
        x2(i+1)=x2(i) + h*f2;
        x2_dot(i) = f2;

    end
    
    % Plot
    subplot(3,1,1)
    plot(t, x1, 'DisplayName', ['h = ' num2str(h)],'LineWidth',2)
    hold on; title("Deslocamento"); ylabel('y (m)'); legend

    subplot(3,1,2)
    plot(t, x2, 'DisplayName', ['h = ' num2str(h)],'LineWidth',2)
    hold on; title("Velocidade"); ylabel('dy/dt (m/s)'); legend

    subplot(3,1,3)
    plot(t, x2_dot, 'DisplayName', ['h = ' num2str(h)],'LineWidth',2)
    hold on; title("Aceleração"); ylabel('d²y/dt² (m/s²)'); xlabel('t (s)'); legend

    sgtitle('Método de Euler')

end



function degrau=z(t)
    if t < 2
        degrau = 0;
    else
        degrau = 0.25;
    end
end