clc, clear

% Criacao da malha

h = 3;
d = 5*h;
L = 2*h;
H = 8*h;

dx = 1;
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


% Plota apenas os pontos da malha (isSolido == false)
mask = ~isSolido;  % true para pontos de malha válidos
scatter( X(mask), Y(mask), 20, 'k', 'filled' )
hold on
scatter(x_galpao, y_galpao, 10, 'r', 'filled')
axis equal
xlabel('x'); ylabel('y');
title('Pontos da malha (excluindo o interior do galpão)');
% Configurar o grid com espaçamento dx e dy
xticks(min(x):dx:max(x));  % Linhas verticais a cada dx
yticks(min(y):dy:max(y));  % Linhas horizontais a cada dy
grid on;

psi = 500*ones(size(X));
psi(isSolido) = 0;
psi = cond_contorno(X,Y,dx,dy,H,V,psi);

erro = 1;
while max(abs(erro(:))) > toleracia
    psi_aux = psi;
    for i=2:size(X,2)-1
        for j=2:size(Y,1)-1
            if ~isSolido(j,i)
                if isSolido(j,i-1) || isSolido(j,i+1) || isSolido(j-1,i)
                    if isSolido(j, i+1) && ~isSolido(j-1, i)
                        %Quadrante 1
                        a = (d-X(j,i))/dx;
                        psi(j,i)=(2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                    
                    elseif isSolido(j, i+1) && isSolido(j-1, i)
                        %Quadrante 2
                        a = ((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx;
                        b = (Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy;
                        psi(j,i)=(2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i-1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1))); 
                    
                    elseif isSolido(j, i-1) && isSolido(j-1, i)
                        %Quadrante 3
                        a = ((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx;
                        b = (Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy;
                        psi(j,i)=(2/(a*dx^2)+2/(b*dy^2))^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+2*psi(j+1,i)/(dy^2*(b+1)));
                        
                    else %isSolido(j, i-1) && ~isSolido(j-1, i)
                        %Quadrante 4
                        a = (X(j,i)-(d+L))/dx;
                        psi(j,i)=(2/(a*dx^2)+2/dy^2)^(-1)*(2*psi(j,i+1)/(dx^2*(a+1))+(psi(j+1,i)+psi(j-1,i))/dy^2);
                        
                    end
                else
                    %Pontos interiores
                    psi(j,i) = (psi(j,i+1)+psi(j,i-1)+psi(j+1,i)+psi(j-1,i))/4;
                end    
            end 
        end
    end
    psi(2:end-1,2:end-1)=lambda*psi(2:end-1,2:end-1)+(1-lambda)*psi_aux(2:end-1,2:end-1);
    psi = cond_contorno(X,Y,dx,dy,H,V,psi);
    erro = abs((psi-psi_aux)./psi);
    display(max(erro(:)));
    contourf(X, Y, psi, 20, 'LineColor', 'none');
    colorbar;
    axis equal;
    xlabel('x'); ylabel('y');
    title('Função de Corrente Genérica (\psi)');
end



% % Exibir a matriz psi - Item A
% figure;
% contourf(X, Y, psi, 20, 'LineColor', 'none');
% colorbar;
% axis equal;
% xlabel('x'); ylabel('y');
% title('Função de Corrente Genérica (\psi)');
% hold on;
% scatter(X(isSolido), Y(isSolido), 2, 'r', 'filled'); % Destacar o sólido
% 
% % Calcular velocidades (apenas para pontos fora do sólido)
% u = zeros(size(psi));
% v = zeros(size(psi));
% u(~isSolido) = gradient(psi(~isSolido), dy);  % u = ∂ψ/∂y
% v(~isSolido) = -gradient(psi(~isSolido), dx); % v = -∂ψ/∂x
% 
% % Plotar campo de velocidade (setas pretas) - ItemB
% figure;
% quiver(X(mask), Y(mask), u(mask), v(mask), ...
%     'AutoScaleFactor', 1.5, ...
%     'Color', 'k', ...          % Define cor preta
%     'LineWidth', 1.2);         % Espessura das setas
% axis equal;
% grid on;                       % Adiciona grade
% xlabel('x (m)'); 
% ylabel('y (m)');
% title('Campo de Velocidade (Genérico)');
% 
% % Ajustar limites dos eixos (se necessário)
% xlim([min(x), max(x)]);
% ylim([min(y), max(y)]);
% 
% %Pressão item C
% 
% % Calcular a pressão (equação corrigida)
% gamma_term = (gamma_ar - 1) / gamma_ar;
% DeltaP = rho * gamma_term * ( (V^2)/2 - (u.^2 + v.^2)/2 );
% 
% % Mascarar pontos dentro do sólido (isSolido)
% DeltaP(isSolido) = NaN; % Ignorar pressão no interior
% 
% % Plotar pressão
% figure;
% contourf(X, Y, DeltaP, 50, 'LineColor', 'none');
% colorbar;
% colormap('jet');
% hold on;
% scatter(X(isSolido), Y(isSolido), 2, 'k', 'filled'); % Destacar o sólido
% axis equal;
% xlabel('x (m)');
% ylabel('y (m)');
% title('Distribuição de Pressão Relativa (p - p_{atm}) [Pa]');
% 
% % Item D

function psi = cond_contorno(X,Y,dx,dy,H,V,psi)
    
    %Margem inferior
    psi(1,:)=0;
    
    %Margem esquerda
    for j = 2:size(Y,1)-1
        if Y(j,1)<=0.05*H
            psi(j,1)=psi(j,2)+V*dx^2/(0.1*H);
        else
            psi(j,1)=psi(j,2);
        end
    end
    
    %Margem direita
    for j = 2:size(Y,1)-1
        if Y(j,end)<=0.05*H
            psi(j,end)=psi(j,end-1)+V*dx^2/(0.1*H);
        else
            psi(j,end)=psi(j,end-1);
        end
    end
    
    %Margem superior
    for i = 1:size(X,2)
        psi(end,i) = psi(end-1,i)+dy*V;
    end

end



