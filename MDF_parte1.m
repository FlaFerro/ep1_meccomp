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

% Plota apenas os pontos da malha (isSolido == false)
mask = ~isSolido;  % true para pontos de malha válidos
scatter( X(mask), Y(mask), 1, 'k', 'filled' ) 
axis equal
xlabel('x'); ylabel('y');
title('Pontos da malha (excluindo o interior do galpão)');

% Resolucao da eq de Laplace
phi = randi([4, 6], length(y), length(x));
erro = 1;
cont = 0;
%while max(abs(erro(:))) > toleracia
while cont < 1
    phi = cond_contorno(x,y,H,V,dx,dy,phi);
    phi_aux = phi;
    for i= 2:length(x)/2+1
        for j = 2:(length(y)-1)
            phi_aux(i,j) = passo(phi,lambda,i,j);
            phi_aux(i,length(x)-j+1) = phi_aux(i,j);
        end
    end
    erro = phi-phi_aux;
    disp(max(abs(erro(:))));
    phi = phi_aux;
    cont = cont + 1;
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
%Aproveitando a simetria
        phi(j,length(x)) = phi(j,1);
    end
%Margem superior e inferior
    for i = 1:length(x)/2 + 1
        %Superior
        phi(length(y),i) = phi(length(y)-1,i) + dy*V;
        %Inferior
        phi(1,i) = 0;
%Aproveitando a simetria
        %Superior
        phi(length(y),length(x)-i+1) =  phi(length(y)-1,i) + dy*V;
        %Inferior
        phi(1,length(x)-i+1) = 0 ;
    end
%Margem inferior
end

function val_phi_novo = passo(phi, lambda,x,y)
    val_phi_novo = (1-lambda)*phi(x,y) + (lambda/4)*(phi(x-1,y)+phi(x+1,y)+phi(x,y-1)+phi(x,y+1));
end