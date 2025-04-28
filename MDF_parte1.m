clc, clear

% Criacao da malha

h = 3;
d = 5*h;
L = 2*h;
H = 8*h;

dx = 0.1;
dy = dx;

x = 0:dx:(2*d+L);
y = 0:dy:H;

[X,Y] = meshgrid(x,y);

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
scatter( X(mask), Y(mask), 1, 'k', 'filled' ) 
axis equal
xlabel('x'); ylabel('y');
title('Pontos da malha (excluindo o interior do galpão)');

% Resolucao da eq de Laplace

phi = zeros(size(X));

%Condicoes de contorno



