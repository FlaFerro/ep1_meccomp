clc, clear

% Criacao da malha

h = 3;
d = 5*h;
L = 2*h;
H = 8*h;

dx = h/8;
dy = dx;

x = 0:dx:(2*d+L);
y = 0:dy:H;

[X,Y] = meshgrid(x,y);

% Outros parâmetros do problema

V = 100*3.6; %Velocidade em m/s!
lambda = 1.85;
toleracia = 0.01;
rho = 1.25;    % Densidade do ar (kg/m³)
gamma_ar = 1.4;% Razão de calor específico do ar

% Verificacao dos pontos dentro do galpao

isSolido = false(size(X));

for i = 1:length(x)
    for j = 1:length(y)
        xi = x(i);
        yj = y(j);
        % Se o ponto está dentro do telhado
        if xi >= d && xi <= d + L
            y_teto = sqrt((L/2)^2 - (xi - (d + L/2))^2) + h;
            if yj <= y_teto
                isSolido(j,i) = true; 
            end
        end
    end
end



% Plota apenas os pontos da malha (isSolido == false)
mask = ~isSolido;  % true para pontos de malha válidos
scatter( X(mask), Y(mask), 2, 'k', 'filled' ) 
axis equal
xlabel('x'); ylabel('y');
title('Pontos da malha (excluindo o interior do galpão)');

% Criação de matriz de testes

% Inicializar psi com escoamento uniforme (V*y)
psi = V * Y; % ψ = V*y (escoamento uniforme na direção x)

% Ajustar condições de contorno:
% 1. ψ = 0 na base (y = 0)
psi(Y == 0) = 0;

% 2. Dentro do sólido, ψ = constante (aqui, 0)
psi(isSolido) = 0;

% 3. Condições laterais (opcional, se necessário):
% - Na entrada (x = 0) e saída (x = 2d + L), manter ψ = V*y
% - Se houver paredes laterais sólidas, ajustar ψ para constante

% Exibir a matriz psi - Item A
figure;
contourf(X, Y, psi, 20, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x'); ylabel('y');
title('Função de Corrente Genérica (\psi)');
hold on;
scatter(X(isSolido), Y(isSolido), 2, 'r', 'filled'); % Destacar o sólido

% Calcular velocidades (apenas para pontos fora do sólido)
u = zeros(size(psi));
v = zeros(size(psi));
u(~isSolido) = gradient(psi(~isSolido), dy);  % u = ∂ψ/∂y
v(~isSolido) = -gradient(psi(~isSolido), dx); % v = -∂ψ/∂x

% Plotar campo de velocidade (setas pretas) - ItemB
figure;
quiver(X(mask), Y(mask), u(mask), v(mask), ...
    'AutoScaleFactor', 1.5, ...
    'Color', 'k', ...          % Define cor preta
    'LineWidth', 1.2);         % Espessura das setas
axis equal;
grid on;                       % Adiciona grade
xlabel('x (m)'); 
ylabel('y (m)');
title('Campo de Velocidade (Genérico)');

% Ajustar limites dos eixos (se necessário)
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

%Pressão item C

% Calcular a pressão (equação corrigida)
gamma_term = (gamma_ar - 1) / gamma_ar;
DeltaP = rho * gamma_term * ( (V^2)/2 - (u.^2 + v.^2)/2 );

% Mascarar pontos dentro do sólido (isSolido)
DeltaP(isSolido) = NaN; % Ignorar pressão no interior

% Plotar pressão
figure;
contourf(X, Y, DeltaP, 50, 'LineColor', 'none');
colorbar;
colormap('jet');
hold on;
scatter(X(isSolido), Y(isSolido), 2, 'k', 'filled'); % Destacar o sólido
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('Distribuição de Pressão Relativa (p - p_{atm}) [Pa]');

% Item D



