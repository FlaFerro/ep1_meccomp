clc, clear

% Criacao da malha

h = 3;
d = 5*h;
L = 2*h;
H = 8*h;

dx = 2;
dy = dx;

x = 0:dx:(2*d+L);
y = 0:dy:H;

[X,Y] = meshgrid(x,y);

% Outros parâmetros do problema

V = 100*3.6; %Velocidade em m/s!
lambda = 1.85;
toleracia = 0.01;

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
phi = zeros(size(X));
erro = 1;
cont = 0;
while cont < 100 %max(abs(erro(:))) > toleracia
    phi = cond_contorno(x,y,H,V,dx,dy,phi);
    phi_aux = phi;
    for coluna= 2:length(x)/2+1 
        for linha = 2:(length(y)-1) 
            phi_aux(linha,coluna) = passo(phi,lambda,linha,coluna,isSolido,d,dx,L,h,dy);
            phi_aux(linha,length(x)-coluna+1) = phi_aux(linha,coluna); %Simetria
        end
    end
    erro = phi-phi_aux;
    disp(max(abs(erro(:))));
    phi = phi_aux;
    cont = cont+1;
    mesh(phi);
end

%Condicoes de contorno

function phi = cond_contorno(x,y,H,V,dx,dy,phi)
%Margem esquerda
    for j = 1:length(y)
        if j*dy <= 0.05*H
             phi(j,1) = V*(dx^2)/(0.1*H) + phi(j,2);
        else
            phi(j,1) = phi(j,2);
        end
        phi(j,length(x)) = phi(j,1); %Aproveitando a simetria
    end

%Margem superior e inferior
    for i = 1:length(x)
        %Superior
        phi(length(y),i) = phi(length(y)-1,i) + dy*V;
        %Inferior
        phi(1,i) = 0;
% %Aproveitando a simetria
            %Superior
       phi(length(y),length(x)-i+1) =  phi(length(y),i);
       %Inferior
       phi(1,length(x)-i+1) = 0 ;
    end
end

function val_phi_novo = passo(phi, lambda,linha,coluna,isSolido,d,dx,L,h,dy)
    
    %Define o valor no prédio.
    if isSolido(linha,coluna) 
        val_phi_novo = 0;  

    %Utiliza cond. de contorno irregulares se perto do prédio.
    elseif isSolido(linha,coluna+1) || isSolido(linha-1,coluna) 
        a = 1;
        b = 1;
        if (coluna-1)*dx < d %Só na parede lateral.
            a = (d-(coluna-1)*dx)/dx ;
        else %Nos contornos do teto circular
            a = ((d+L/2-(coluna-1)*dx) - sqrt( (L/2)^2 - ((linha-1)*dy - h)^2 ))/dy;
            if ~isSolido(linha,coluna+1)
                a=1;
            end
            b = ((linha-1)*dy -h - sqrt((L/2)^2 - (d+(L/2) - (coluna-1)*dx)^2))/dy;
            if ~isSolido(linha-1,coluna)
                b=1;
            end
        end
        %Calcula o valor usando as proporções.
        val_phi_novo = ((a*b*dx^2*dy^2)/(a*dx^2 + b*dy^2))*((phi(linha,coluna-1)/(dx^2*(a+1)))+(phi(linha+1,coluna)/dy^2*(b+1)));
    else
    %Executa sobrerrelaxação tradicional se não está perto do prédio.
    val_phi_novo = (1-lambda)*phi(linha,coluna) + (lambda/4)*(phi(linha-1,coluna)+phi(linha+1,coluna)+phi(linha,coluna-1)+phi(linha,coluna+1));
    end
end