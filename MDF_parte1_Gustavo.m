clc, clear

% Criacao da malha

h = 3;
d = 5*h;
L = 2*h;
H = 8*h;

dx = 0.15;
dy = dx;

x = 0:dx:(2*d+L);
y = 0:dy:H;

[X,Y] = meshgrid(x,y);

% Outros parâmetros do problema

V = 100/3.6; %Velocidade em m/s!
lambda = 1.85;
toleracia = 0.001;
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

[~, j_first] = find(isSolido(1, :), 1, 'first');
x_telha = (X(1,j_first)):dx:(d + L - dx);
y_telha = sqrt((L/2).^2 - (x_telha - (d + L/2)).^2) + h;
y_telha_aux = Y(length(0:dy:h)+1,1):dy:(h+L/2);
x_telha_esq = d + L/2 - sqrt((L/2)^2 - (y_telha_aux - h).^2); % Ramo esquerdo (sinal -)
x_telha_dir  = d + L/2 + sqrt((L/2)^2 - (y_telha_aux - h).^2); % Ramo direito (sinal +)
y_parede_esq = 0:dy:h;
x_parede_esq = d*ones(size(y_parede_esq));
y_parede_dir = 0:dy:h;
x_parede_dir = (d+L)*ones(size(y_parede_dir));

x_galpao = [x_parede_esq, x_telha_esq, x_telha, x_telha_dir, x_parede_dir];
y_galpao = [y_parede_esq, y_telha_aux, y_telha, y_telha_aux, y_parede_dir];


% % Plota apenas os pontos da malha (isSolido == false)
mask = ~isSolido;  % true para pontos de malha válidos
% scatter( X(mask), Y(mask), 20, 'k', 'filled' )
% hold on
% scatter(x_galpao, y_galpao, 10, 'r', 'filled')
% axis equal
% xlabel('x'); ylabel('y');
% title('Pontos da malha (excluindo o interior do galpão)');
% % Configurar o grid com espaçamento dx e dy
% xticks(min(x):dx:max(x));  % Linhas verticais a cada dx
% yticks(min(y):dy:max(y));  % Linhas horizontais a cada dy
% grid on;

psi = 0*ones(size(X));
psi(isSolido) = 0;
psi = cond_contorno(X,Y,dx,dy,H,V,psi);
psi_aux=psi;

erro = 1;
while max(abs(erro(:))) > toleracia
    mesh(X,Y,psi)
    psi_aux=psi;
    for i=2:size(X,2)-1
        for j=2:size(Y,1)-1
            if ~isSolido(j,i)
                if isSolido(j,i-1) || isSolido(j,i+1) || isSolido(j-1,i)
                    if isSolido(j, i+1) && ~isSolido(j-1, i)
                        
                        if Y(j,i) > h
                            a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs((d-X(j,i))/dx);
                        end
                        psi(j,i)=(2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                   
                    elseif isSolido(j, i+1) && isSolido(j-1, i)
                        
                        a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i)=(2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1)));
                   
                    elseif isSolido(j, i-1) && isSolido(j-1, i)
                        
                        a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i)=(2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1)));
                       
                    elseif isSolido(j, i-1) && ~isSolido(j-1, i)
                        
                        if Y(j,i) > h
                            a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs(((d+L)-X(j,i))/dx);
                        end
                        psi(j,i)=(2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                       
                    elseif (~isSolido(j,i-1) && ~isSolido(j,i+1)) && isSolido(j-1,i)
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        psi(j,i) = (2/(dx^2) + 2/(b*dy^2))^(-1)*((psi(j,i+1)+psi(j,i-1))/(dx^2)+2*psi(j+1,i)/(b+1));
                    else
                        disp("não devia ter entrado aqui")
                    end
                else
                    %Pontos interiores
                    psi(j,i) = (lambda/4)*(psi(j-1,i)+psi(j+1,i)+psi(j,i-1)+psi(j,i+1)) + (1-lambda)*psi_aux(j,i);
                end    
            end
        end
    end
    erro = abs((psi - psi_aux) ./ max(abs(psi), 1e-10));
    display(max(erro(:)));
end

figure;
contourf(X, Y, psi, 20, 'LineColor', 'none');
cb = colorbar;
ylabel(cb, 'Função de corrente (m^2/s)');
axis equal;
xlabel('x'); ylabel('y');
title('Função de Corrente (\psi)');
hold on;
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled'); % Destacar o sólido

% Calcular velocidades (apenas para pontos fora do sólido)
u = zeros(size(psi));
v = zeros(size(psi));
for linha = 2:length(y)-1
    for coluna = 2:length(x)-1
        v(linha,coluna) = -(psi(linha,coluna+1) - psi(linha,coluna-1)) / (2*dx);
        u(linha,coluna) = (psi(linha+1,coluna) - psi(linha-1,coluna)) / (2*dy);
    end
end

% Calcular pressão no domínio

p = zeros(size(psi));
for linha = 2:length(y)-1
    for coluna = 2:length(x)-1
        p(linha,coluna) = rho * (gamma_ar-1)/gamma_ar * ( (V^2)/2 - ((sqrt(u(linha,coluna)^2 + v(linha,coluna)^2 ))^2)/2) ;
    end
end

% Calcular a pressão ao longo do telhado.
p_telhado = zeros(2, L/dx+1);
for coluna = 1:length(p_telhado)
    p_telhado(1,coluna) = min(p(:,coluna + d/dx - 1));
    p_telhado(2,coluna) = d + (coluna-1)*dx;
end

figure;
yyaxis right
x_telhado_plot = d:dx:(d+L);
plot(x_telhado_plot,sqrt((L/2).^2 - (x_telhado_plot - (d + L/2)).^2) + h , 'Color', 'k', 'LineWidth', 1.5); % Telhado em preto
ylabel('Y do telhado (m)', 'Color', 'k');
ylim([min(y_telha) max(y_telha)]);
set(gca, 'YColor', 'k'); % Cor do eixo esquerdo

yyaxis left
plot(p_telhado(2,:), p_telhado(1,:), 'Color', 'b', 'LineWidth', 1.5); % Pressão em azul escuro
ylabel('Diferencial de pressão (Pa)', 'Color', 'b');
ylim([min(p_telhado(1,:)) max(p_telhado(1,:))]);
set(gca, 'YColor', 'b'); % Cor do eixo direito

xlabel('x (m)');
title('Diferencial de pressão ao longo do telhado');
grid on;

x_c = d + L/2;  % Centro do telhado curvo
y_c = h;
p_forca=p_telhado(:,2:end-1);
i = 1;
F=0;
for xi = d+dx:dx:(d+L)-dx
    yi = sqrt((L/2)^2 - (xi - (d + L/2)).^2) + h;
    n_vec = [xi - x_c, yi - y_c];
    norm_n = norm(n_vec);
    n = n_vec / norm_n;
    
    x_seg=xi+dx;
    y_seg = sqrt((L/2)^2 - (x_seg - (d + L/2)).^2) + h;
    
    dl = sqrt((xi-x_seg)^2+(yi-y_seg)^2);
    dF = -p_forca(1,i)*dl*60*n;
    F = F+dF(2);
    i=i+1;
    
end

fprintf('Força telhado: F = %.2f kN\n', F/1000);

% Plotar campo de velocidade - ItemB
step = 10;
% Cria máscara esparsa
Xq = X(1:step:end, 1:step:end);
Yq = Y(1:step:end, 1:step:end);
uq = u(1:step:end, 1:step:end);
vq = v(1:step:end, 1:step:end);
maskq=mask(1:step:end, 1:step:end);
vel_mag = sqrt(u.^2 + v.^2);
figure;
contourf(X(2:end-1,2:end-1), Y(2:end-1,2:end-1), vel_mag(2:end-1,2:end-1), 20, 'LineColor', 'none'); % 20 níveis de contorno
cb = colorbar;
ylabel(cb, 'Velocidade (m/s)');  % Adiciona legenda ao colorbar
colormap jet; % Mapa de cores (você pode trocar, por ex: 'parula', 'hot', etc)
hold on;
quiver(Xq(maskq), Yq(maskq), uq(maskq), vq(maskq), ...
    'AutoScaleFactor', 1.2, ...
    'Color', 'k', ...          % Define cor preta
    'LineWidth', 1);         % Espessura das setas
axis equal;
grid on;                       % Adiciona grade
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled'); % Destacar o sólido
xlabel('x (m)');
ylabel('y (m)');
title('Campo de Velocidade');


figure;
contourf(X, Y, p, 12 , 'LineColor', 'none');
cb = colorbar;
ylabel(cb, 'Pressão (Pa)');
axis equal;
xlabel('x'); ylabel('y');
title('Distribuição diferencial de pressão');
hold on; 
scatter(X(isSolido), Y(isSolido), 10, 'k', 'filled'); % Destacar o sólido

function psi = cond_contorno(X,Y,dx,dy,H,V,psi)
   
    %Margem inferior
    psi(1,:)=0;
   
    %Margem esquerda e direita
    for j = 2:size(Y,1)
        if Y(j,1)<=0.05*H
            psi(j,1)=V*(Y(j,1))^2/(0.1*H);
            psi(j,end)=V*(Y(j,end))^2/(0.1*H);
        else
            psi(j,1)=V*Y(j,1)+V*(1-0.05*H);
            psi(j,end)=V*Y(j,end)+V*(1-0.05*H);
        end
    end
   
    psi(end,:) = psi(end,1);

end