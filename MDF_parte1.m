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

V = 100/3.6; %Velocidade em m/s!
lambda = 1.85;
toleracia = 0.1;

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

% Resolucao da eq de Laplace
psi = zeros(size(X));
psi_aux = psi;
erro = 1;
psi = cond_contorno(x,y,H,V,dx,dy,psi);
while max(abs(erro(:))) > toleracia
    for coluna= 2:length(x)/2+1 
        for linha = 2:(length(y)-1) 
            psi(linha,coluna) = passo(psi,lambda,linha,coluna,isSolido,d,dx,L,h,dy);
            psi(linha,length(x)-coluna+1) = psi(linha,coluna); %Simetria
        end
    end
    erro = psi-psi_aux;
    disp(max(abs(erro(:))));
    psi_aux = psi;
end


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


%Condicoes de contorno

function psi = cond_contorno(x,y,H,V,dx,dy,psi)

% %Aproveitando a simetria
            %Superior
%        psi(length(y),length(x)-i+1) =  psi(length(y),i);
       %Inferior
    %   psi(1,length(x)-i+1) = 0 ;
    %end

    %Margem esquerda
    for j = 2:length(y)
        psi(j,1) = V*(j-1)*dx - 0.0025*V*H;
        psi(j,length(x)) = psi(j,1); %Aproveitando a simetria
    end
      for i = 2:length(x)
        psi(length(y),i) = psi(length(y),1);
      end

        %Inferior
        psi(1,i) = 0;

 end

function val_psi_novo = passo(psi, lambda,linha,coluna,isSolido,d,dx,L,h,dy)
    
    %Define o valor no prédio.
    if isSolido(linha,coluna) 
        val_psi_novo = 0;  

    %Utiliza cond. de contorno irregulares se perto do prédio.
    elseif isSolido(linha,coluna+1) || isSolido(linha-1,coluna) 
        a = 1;
        b = 1;
        if (coluna-1)*dx < d %Só na parede lateral.
            a = (d-(coluna-1)*dx)/dx ;
        else %Nos contornos do teto circular
            a = ((d + L/2 -(coluna-1)*dx) - sqrt( (L/2)^2 - ((linha-1)*dy - h)^2 ))/dx;
            if ~isSolido(linha,coluna+1)
                a=1;
            end
            b = ((linha-1)*dy -h - sqrt((L/2)^2 - (d+(L/2) - (coluna-1)*dx)^2))/dy;
            if ~isSolido(linha-1,coluna)
                b=1;
            end
        end
        %Calcula o valor usando as proporções.
        val_psi_novo = ((a*b*dx^2*dy^2)/(a*dx^2 + b*dy^2))*((psi(linha,coluna-1)/(dx^2*(a+1)))+(psi(linha+1,coluna)/dy^2*(b+1)));
    else
    %Executa sobrerrelaxação tradicional se não está perto do prédio.
    val_psi_novo = (lambda/4)*(psi(linha-1,coluna)+psi(linha+1,coluna)+psi(linha,coluna-1)+psi(linha,coluna+1)) + (1-lambda)*psi(linha,coluna);
    end
end