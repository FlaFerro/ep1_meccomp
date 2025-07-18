clc, clear

% Criacao da malha

h = 3;          % Altura da parede do galpao
d = 5*h;        % Distancia da origem ate o galpao
L = 2*h;        % Largura do galpao
H = 8*h;        % Altura total do dominio

dx = 0.5;       % Passo em x
dy = dx;        % Passo em y (igual a dx)

x = 0:dx:(2*d+L);   % Vetor de posicoes em x
y = 0:dy:H;         % Vetor de posicoes em y

[X,Y] = meshgrid(x,y);  % Criar malha 2D

% Parametros do problema

V = 100/3.6;        % Velocidade do escoamento (m/s)
lambda = 1.85;      % Parametro de relaxacao
tolerancia = 0.01;  % Tolerancia para convergencia
rho = 1.25;         % Densidade do ar (kg/m³)
gamma_ar = 1.4;     % Razao de calor especifico do ar

% Identificar pontos dentro do galpao (solido)

isSolido = false(size(X));  % Matriz logica para identificar regiao solida

for i = 1:length(x)
    for j = 1:length(y)
        xi = x(i);
        yj = y(j);
        % Verifica se o ponto esta dentro do telhado curvo
        if xi >= d && xi <= d + L
            y_teto = sqrt((L/2)^2 - (xi - (d + L/2))^2) + h;
            if yj <= y_teto
                isSolido(j,i) = true;
            end
        end
    end
end

% Definir geometria do galpao para plotagem

[~, j_first] = find(isSolido(1, :), 1, 'first');
x_telha = (X(1,j_first)):dx:(d + L - dx);
y_telha = sqrt((L/2).^2 - (x_telha - (d + L/2)).^2) + h;
y_telha_aux = Y(length(0:dy:h)+1,1):dy:(h+L/2);
x_telha_esq = d + L/2 - sqrt((L/2)^2 - (y_telha_aux - h).^2);
x_telha_dir  = d + L/2 + sqrt((L/2)^2 - (y_telha_aux - h).^2);
y_parede_esq = 0:dy:h;
x_parede_esq = d*ones(size(y_parede_esq));
y_parede_dir = 0:dy:h;
x_parede_dir = (d+L)*ones(size(y_parede_dir));

x_galpao = [x_parede_esq, x_telha_esq, x_telha, x_telha_dir, x_parede_dir];
y_galpao = [y_parede_esq, y_telha_aux, y_telha, y_telha_aux, y_parede_dir];

% Plotar malha computacional
figure;
mask = ~isSolido;  % Pontos fluidos
scatter( X(mask), Y(mask), 20, 'k', 'filled' )
hold on
scatter(x_galpao, y_galpao, 10, 'r', 'filled')
axis equal
xlabel('x'); ylabel('y');
title('Pontos da malha (excluindo o interior do galpao)');
xticks(min(x):dx:max(x));
yticks(min(y):dy:max(y));
grid on;

% Inicializar funcao de corrente

psi = 0*ones(size(X));
psi(isSolido) = 0;
psi = cond_contorno(X,Y,dx,dy,H,V,psi);
psi_aux = psi;

% Resolver equacao de Laplace para psi com metodo SOR

erro = 1;
while max(abs(erro(:))) > tolerancia
%     mesh(X,Y,psi)
    psi_aux = psi;
    
    for i = 2:size(X,2)-1
        for j = 2:size(Y,1)-1
            if ~isSolido(j,i)
                % Tratar pontos proximos as fronteiras
                if isSolido(j,i-1) || isSolido(j,i+1) || isSolido(j-1,i)
                    % Casos especiais para fronteiras curvas
                    if isSolido(j, i+1) && ~isSolido(j-1, i)
                        if Y(j,i) > h
                            a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs((d-X(j,i))/dx);
                        end
                        psi(j,i) = (2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                    
                    elseif isSolido(j, i+1) && isSolido(j-1, i)
                        a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i) = (2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1)));
                    
                    elseif isSolido(j, i-1) && isSolido(j-1, i)
                        a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i) = (2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1)));
                       
                    elseif isSolido(j, i-1) && ~isSolido(j-1, i)
                        if Y(j,i) > h
                            a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs(((d+L)-X(j,i))/dx);
                        end
                        psi(j,i) = (2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                       
                    elseif (~isSolido(j,i-1) && ~isSolido(j,i+1)) && isSolido(j-1,i)
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i) = (2/(dx^2) + 2/(b*dy^2))^(-1)*((psi(j,i+1)+psi(j,i-1))/(dx^2)+2*psi(j+1,i)/(dy^2*(b+1)));
                    else
                        disp("Ponto nao esperado")
                    end
                else
                    % Pontos internos - equacao de Laplace
                    psi(j,i) = (lambda/4)*(psi(j-1,i)+psi(j+1,i)+psi(j,i-1)+psi(j,i+1)) + (1-lambda)*psi_aux(j,i);
                end    
            end
        end
    end
    erro = abs((psi - psi_aux) ./ max(abs(psi), 1e-10));
    display(max(erro(:)));
end

% Plotar funcao de corrente

figure;
contourf(X, Y, psi, 20, 'LineColor', 'none');
cb = colorbar;
ylabel(cb, 'Funcao de corrente (m^2/s)');
axis equal;
xlabel('x'); ylabel('y');
title('Funcao de Corrente (\psi)');
hold on;
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled');

% Calcular campo de velocidade

u = zeros(size(psi));
v = zeros(size(psi));
for linha = 2:length(y)-1
    for coluna = 2:length(x)-1
        v(linha,coluna) = -(psi(linha,coluna+1) - psi(linha,coluna-1)) / (2*dx);
        u(linha,coluna) = (psi(linha+1,coluna) - psi(linha-1,coluna)) / (2*dy);
    end
end
u(length(u(:,1)),:) = V;
for coluna = 1:length(X)
    u(1,coluna) = (psi(2,coluna)-psi(1,coluna))/dy;
end
u(1:end,1) = V;
u(1:end,length(u(1,:))) = V;



% Calcular campo de pressao

p = zeros(size(psi));
for linha = 2:length(y)-1
    for coluna = 2:length(x)-1
        p(linha,coluna) = rho * ((gamma_ar-1)/gamma_ar) * ( (V^2)/2 - ((sqrt(u(linha,coluna)^2 + v(linha,coluna)^2 ))^2)/2 );
    end
end

% Calcular pressao ao longo do telhado

p_telhado = zeros(2, L/dx+1);
for coluna = 1:length(p_telhado)
    p_telhado(1,coluna) = min(p(:,coluna+length(0:dx:d)-1));
    p_telhado(2,coluna) = d + (coluna-1)*dx;
end

% Plotar pressao no telhado

figure;
yyaxis right
x_telhado_plot = d:dx:(d+L);
plot(x_telhado_plot,sqrt((L/2).^2 - (x_telhado_plot - (d + L/2)).^2) + h , 'Color', 'k', 'LineWidth', 1.5);
ylabel('Y do telhado (m)', 'Color', 'k');
ylim([min(y_telha) max(y_telha)]);
set(gca, 'YColor', 'k');

yyaxis left
plot(p_telhado(2,:), p_telhado(1,:), 'Color', 'b', 'LineWidth', 1.5);
ylabel('Diferencial de pressao (Pa)', 'Color', 'b');
ylim([min(p_telhado(1,:)) max(p_telhado(1,:))]);
set(gca, 'YColor', 'b');

xlabel('x (m)');
title('Diferencial de pressao ao longo do telhado');
grid on;

% Calcular forca resultante no telhado

x_c = d + L/2;  % Centro do telhado
y_c = h;
p_forca = p_telhado(:,2:end-1);
i = 1;
F = 0;
for xi = d+dx:dx:(d+L)-dx
    yi = sqrt((L/2)^2 - (xi - (d + L/2)).^2) + h;
    n_vec = [xi - x_c, yi - y_c];
    norm_n = norm(n_vec);
    n = n_vec / norm_n;
    
    x_seg = xi+dx;
    y_seg = sqrt((L/2)^2 - (x_seg - (d + L/2)).^2) + h;
    
    dl = sqrt((xi-x_seg)^2+(yi-y_seg)^2);
    dF = -p_forca(1,i)*dl*60*n;
    F = F+dF(2);
    i = i+1;
end

fprintf('Forca telhado: F = %.2f kN\n', F/1000);

% Plotar campo de velocidade

step = 2;
Xq = X(1:step:end, 1:step:end);
Yq = Y(1:step:end, 1:step:end);
uq = u(1:step:end, 1:step:end);
vq = v(1:step:end, 1:step:end);
maskq = mask(1:step:end, 1:step:end);
vel_mag = sqrt(u.^2 + v.^2);

figure;
contourf(X, Y, vel_mag, 20, 'LineColor', 'none');
cb = colorbar;
ylabel(cb, 'Velocidade (m/s)');
colormap jet;
hold on;
quiver(Xq(maskq), Yq(maskq), uq(maskq), vq(maskq), ...
    'AutoScaleFactor', 1.2, ...
    'Color', 'k', ...
    'LineWidth', 1);
axis equal;
grid on;
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled');
xlabel('x (m)');
ylabel('y (m)');
title('Campo de Velocidade');

% Plotar campo de pressao

figure;
contourf(X(2:end-1,2:end-1), Y(2:end-1,2:end-1), p(2:end-1,2:end-1), 12 , 'LineColor', 'none');
cb = colorbar;
ylabel(cb, 'Pressao (Pa)');
axis equal;
xlabel('x'); ylabel('y');
title('Distribuicao diferencial de pressao');
hold on; 
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled');

% Funcao para condicoes de contorno

function psi = cond_contorno(X,Y,dx,dy,H,V,psi)
    % Margem inferior
    psi(1,:) = 0;
   
    % Margens laterais
    for j = 2:size(Y,1)
        if Y(j,1) <= 0.05*H
            psi(j,1) = V*(Y(j,1))^2/(0.1*H);
            psi(j,end) = V*(Y(j,end))^2/(0.1*H);
        else
            psi(j,1) = V*Y(j,1) + V*(1-0.05*H);
            psi(j,end) = V*Y(j,end) + V*(1-0.05*H);
        end
    end
   
    psi(end,:) = psi(end,1);
end