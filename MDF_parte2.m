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

V = 100/3.6; %Velocidade em m/s!
lambda = 1.15;
toleracia = 0.01;
rho = 1.25;    % Densidade do ar (kg/m³)
gamma_ar = 1.4;% Razão de calor específico do ar
T_f = 20;
T_d = 40;
k = 0.026;
c = 1002;

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

load('dados_campo.mat');  % Isso carrega psi, u e v no workspace

T = 20*ones(size(psi));
T(isSolido) = T_d;
T = cond_contorno_temp(X,Y, u, v, T, isSolido, k, dx, dy, rho, c, T_f, T_d, V, d, L);

erro = 1;
while max(abs(erro(:))) > toleracia
    T_aux = T;
    for i = 2:size(X,2)-1
        for j = 2:size(Y,1)-1
            if ~isSolido(j,i)
                if isSolido(j,i-1) || isSolido(j,i+1) || isSolido(j-1,i)
                    if isSolido(j, i+1) && ~isSolido(j-1, i)
                        %Quadrante 1
                        if Y(j,i) > h
                            a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs((d-X(j,i))/dx);
                        end
                        if u(j,i)>=0 && v(j,i)>=0
                            T(j,i)=(rho*c*u(j,i)/dx+rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(u(j,i)*T(j,i-1)/dx+v(j,i)*T(j-1,i)/dy));
                        elseif u(j,i)>=0 && v(j,i)<0
                            T(j,i)=(rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(u(j,i)*T(j,i-1)/dx-v(j,i)*T(j+1,i)/dy));
                        elseif u(j,i)<0 && v(j,i)>=0
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)+rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(-u(j,i)*T_d/(a*dx)+v(j,i)*T(j-1,i)/dy));
                        else
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(-u(j,i)*T_d/(a*dx)-v(j,i)*T(j+1,i)/dy));
                        end
                        
                    elseif isSolido(j, i+1) && isSolido(j-1, i)
                        %Quadrante 2
                        a = abs(((d + L/2 - sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        if u(j,i)>=0 && v(j,i)>=0
                            T(j,i)=(rho*c*u(j,i)/dx+rho*c*v(j,i)/(b*dy)+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx+v(j,i)*T_d/(b*dy)));
                        elseif u(j,i)>=0 && v(j,i)<0
                            T(j,i)=(rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx-v(j,i)*T(j+1,i)/dy));
                        elseif u(j,i)<0 && v(j,i)>=0
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)+rho*c*v(j,i)/(b*dy)+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(-u(j,i)*T_d/(a*dx)+v(j,i)*T_d/(b*dy)));
                        else
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i-1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T_d/(a*dx)-v(j,i)*T(j+1,i)/dy));
                        end
                        
                    elseif isSolido(j, i-1) && isSolido(j-1, i)
                        %Quadrante 3
                        a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        if u(j,i)>=0 && v(j,i)>=0
                            T(j,i)=(rho*c*u(j,i)/dx+rho*c*v(j,i)/(b*dy)+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx+v(j,i)*T_d/(b*dy)));
                        elseif u(j,i)>=0 && v(j,i)<0
                            T(j,i)=(rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx-v(j,i)*T(j+1,i)/dy));
                        elseif u(j,i)<0 && v(j,i)>=0
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)+rho*c*v(j,i)/(b*dy)+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(-u(j,i)*T_d/(a*dx)+v(j,i)*T_d/(b*dy)));
                        else
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(b*dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+2/(dy^2)*(T_d/(b*(b+1))+T(j+1,i)/(b+1)))+rho*c*(u(j,i)*T_d/(a*dx)-v(j,i)*T(j+1,i)/dy));
                        end
                              
                    elseif isSolido(j, i-1) && ~isSolido(j-1, i)
                        %Quadrante 4
                        if Y(j,i) > h
                            a = abs(((d + L/2 + sqrt((L/2)^2 - (Y(j,i) - h)^2))-X(j,i))/dx);
                        else
                            a = abs(((d+L)-X(j,i))/dx);
                        end
                        
                        if u(j,i)>=0 && v(j,i)>=0
                            T(j,i)=(rho*c*u(j,i)/dx+rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(u(j,i)*T(j,i-1)/dx+v(j,i)*T(j-1,i)/dy));
                        elseif u(j,i)>=0 && v(j,i)<0
                            T(j,i)=(rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(u(j,i)*T(j,i-1)/dx-v(j,i)*T(j+1,i)/dy));
                        elseif u(j,i)<0 && v(j,i)>=0
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)+rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(-u(j,i)*T_d/(a*dx)+v(j,i)*T(j-1,i)/dy));
                        else
                            T(j,i)=(-rho*c*u(j,i)/(a*dx)-rho*c*v(j,i)/dy+2*k/(a*dx^2)+2*k/(dy^2))^(-1)*(k*(2/(dx^2)*(T_d/(a*(a+1))+T(j,i+1)/(a+1))+(T(j+1,i)+T(j-1,i))/dy^2)+rho*c*(-u(j,i)*T_d/(a*dx)-v(j,i)*T(j+1,i)/dy));
                        end
                        
                    elseif (~isSolido(j,i-1) && ~isSolido(j,i+1)) && isSolido(j-1,i)
                        b = abs((Y(j,i)-(sqrt((L/2)^2 - (X(j,i) - (d + L/2))^2) + h))/dy);
                        
                        if u(j,i)>=0 && v(j,i)>=0
                            T(j,i)=(2*k/dx^2+k/(b*dy^2)+rho*c*u(j,i)+rho*c*v(j,i)/(b*dy))^(-1)*(k*((T(j,i+1)+T(j,i-1))/dx^2+(2/dy^2)*(T_d/(b*(b+1))+T(j,i-1)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx+v(j,i)*T_d/(b*dy))); 
                        elseif u(j,i)>=0 && v(j,i)<0
                            T(j,i)=(2*k/dx^2+k/(b*dy^2)+rho*c*u(j,i)-rho*c*v(j,i)/(dy))^(-1)*(k*((T(j,i+1)+T(j,i-1))/dx^2+(2/dy^2)*(T_d/(b*(b+1))+T(j,i-1)/(b+1)))+rho*c*(u(j,i)*T(j,i-1)/dx-v(j,i)*T(j+1,i)/dy));
                        elseif u(j,i)<0 && v(j,i)>=0
                            T(j,i)=(2*k/dx^2+k/(b*dy^2)-rho*c*u(j,i)+rho*c*v(j,i)/(b*dy))^(-1)*(k*((T(j,i+1)+T(j,i-1))/dx^2+(2/dy^2)*(T_d/(b*(b+1))+T(j,i-1)/(b+1)))+rho*c*(-u(j,i)*T(j,i+1)/dx+v(j,i)*T_d/(b*dy))); 
                        else
                            T(j,i)=(2*k/dx^2+k/(dy^2)-rho*c*u(j,i)-rho*c*v(j,i)/(b*dy))^(-1)*(k*((T(j,i+1)+T(j,i-1))/dx^2+(2/dy^2)*(T_d/(b*(b+1))+T(j,i-1)/(b+1)))+rho*c*(-u(j,i)*T(j,i+1)/dx-v(j,i)*T(j+1,i)/dy)); 
                        end
                        
                        
                    end
                end
                
                if u(j,i)>=0 && v(j,i)>=0
                    T(j,i)=(rho*c*u(j,i)/dx+rho*c*v(j,i)/dy+2*k/(dx^2)+2*k/(dy^2))^(-1)*(rho*c*u(j,i)*T(j,i-1)/dx+rho*c*v(j,i)*T(j-1,i)/dy+k*(T(j,i+1)+T(j,i-1))/(dx^2)+k*(T(j+1,i)+T(j-1,i))/(dy^2));
                elseif u(j,i)>=0 && v(j,i)<0
                    T(j,i)=(rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(dx^2)+2*k/(dy^2))^(-1)*(rho*c*u(j,i)*T(j,i-1)/dx-rho*c*v(j,i)*T(j+1,i)/dy+k*(T(j,i+1)+T(j,i-1))/(dx^2)+k*(T(j+1,i)+T(j-1,i))/(dy^2));
                elseif u(j,i)<0 && v(j,i)>=0
                    T(j,i)=(-rho*c*u(j,i)/dx+rho*c*v(j,i)/dy+2*k/(dx^2)+2*k/(dy^2))^(-1)*(-rho*c*u(j,i)*T(j,i+1)/dx+rho*c*v(j,i)*T(j-1,i)/dy+k*(T(j,i+1)+T(j,i-1))/(dx^2)+k*(T(j+1,i)+T(j-1,i))/(dy^2));
                else %u(j,i)<0 && v(j,i)<0
                    T(j,i)=(-rho*c*u(j,i)/dx-rho*c*v(j,i)/dy+2*k/(dx^2)+2*k/(dy^2))^(-1)*(-rho*c*u(j,i)*T(j,i+1)/dx-rho*c*v(j,i)*T(j+1,i)/dy+k*(T(j,i+1)+T(j,i-1))/(dx^2)+k*(T(j+1,i)+T(j-1,i))/(dy^2));
                end 
            end     
        end
    end
    T(2:end-1,2:end-1) = lambda*T(2:end-1,2:end-1)+(1-lambda)*T_aux(2:end-1,2:end-1);
    T = cond_contorno_temp(X,Y, u, v, T, isSolido, k, dx, dy, rho, c, T_f, T_d, V, d, L);
    erro = abs((T-T_aux)./T);
    erro = max(erro(:))
end

figure;
contourf(X, Y, T, 20, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x'); ylabel('y');
title('Temperatura');


dTdx = zeros(size(T));
dTdy = zeros(size(T));
for linha = 2:size(T,1)-1
    for coluna = 2:size(T,2)-1
        dTdx(linha,coluna) = (T(linha,coluna+1) - T(linha,coluna-1)) / (2*dx);
        dTdy(linha,coluna) = (T(linha+1,coluna) - T(linha-1,coluna)) / (2*dy);
    end
end

% Qdot = 0;  % Inicializa a taxa de calor total
% x_c = d + L/2;  % Centro do telhado curvo
% y_c=h;
% tol = 1e-6;      % Tolerância para evitar divisão por zero
% 
% for i = 2:(size(X,2)-1)
%     for j = 2:(size(Y,1)-1)
%         if isSolido(j,i-1) || isSolido(j,i+1) || isSolido(j-1,i)
%             x_i = X(j,i); 
%             y_j = Y(j,i);
%             
%             % Determina dl conforme a região
%             if (X(j,i) <= d && Y(j,i) <= h) 
%                 dl = dy;  % Parede vertical
%                 n = [-1, 0]; % Normal para direita ou esquerda (ajuste conforme necessário)
%             elseif (X(j,i) >= d+L && Y(j,i) <= h) 
%                 dl = dy;  % Parede vertical
%                 n = [1, 0]; 
%             else
%                 dl=sqrt(2)*dx;
%                 n_vec = [x_i - x_c, y_j - y_c];
%                 n = n_vec / norm(n_vec);
%             end
%             
%             % Produto escalar gradT · n
%             dotprod = dTdx(j,i)*n(1) + dTdy(j,i)*n(2);
%             
%             % Contribuição de calor
%             Qdot = Qdot - k * dotprod * dl * 60; % 60m de comprimento
%         end
%     end
% end
% 
% fprintf('Taxa de calor perdida: Q = %.2f W\n', Qdot);


figure;
quiver(X(mask), Y(mask), dTdx(mask), dTdy(mask), ...
    'AutoScaleFactor', 1.5, ...
    'Color', 'k', ...          % Define cor preta
    'LineWidth', 1.2);         % Espessura das setas
axis equal;
grid on;                       % Adiciona grade
xlabel('x (m)');
ylabel('y (m)');
title('Gradiente de temperatura');


function T = cond_contorno_temp(X,Y, u, v, T, isSolido, k, dx, dy, rho, c, T_f, T_d, V, d, L)

    %Margem esquerda
    T(:,1) = T_f;
    
    %Margem inferior
    for i = 2:size(X,2)-1
        if ~isSolido(1,i+1) && ~isSolido(1, i-1)
            if u(1,i) >= 0
                T(1,i) = (2*k/dx^2 + 2*k/dy^2 + rho*c*u(1,i))^(-1)*(k*((T(1,i+1)+T(1,i-1))/(dx^2) + 2*T(2,i)/(dy^2)) + rho*c*u(1,i)*T(1,i-1)/dx);                
            else
                T(1,i) = (2*k/dx^2 + 2*k/dy^2 - rho*c*u(1,i))^(-1)*(k*((T(1,i+1)+T(1,i-1))/(dx^2) + 2*T(2,i)/(dy^2)) - rho*c*u(1,i)*T(1,i+1)/dx);
            end
        elseif ((isSolido(1,i+1) && ~isSolido(1, i-1)) || (~isSolido(1,i+1) && isSolido(1, i-1))) && ~isSolido(1,i)
           if isSolido(1,i+1)
               a = (d-X(1,i))/dx;
               if u(1,i) >= 0
                   T(1,i) = (2*k/(a*dx^2)+2*k/dy^2+rho*c*u(1,i)/dx)^(-1)*(k*((2/dx^2)*(T_d/(a*(a+1))+T(1,i-1)/(a+1))+2*T(2,i)/dy^2)+rho*c*u(1,i)*T(1,i-1)/dx);
               else
                   T(1,i) = (2*k/(a*dx^2)+2*k/dy^2-rho*c*u(1,i)/(a*dx))^(-1)*(k*((2/dx^2)*(T_d/(a*(a+1))+T(1,i-1)/(a+1))+2*T(2,i)/dy^2)-rho*c*u(1,i)*T_d/(a*dx));
               end
           else
               a = (X(1,i)-(d+L))/dx;
               if u(1,i) >= 0
                   T(1,i) = (2*k/(a*dx^2)+2*k/dy^2+rho*c*u(1,i)/dx)^(-1)*(k*((2/dx^2)*(T_d/(a*(a+1))+T(1,i+1)/(a+1))+2*T(2,i)/dy^2)+rho*c*u(1,i)*T(1,i-1)/dx);
               else
                   T(1,i) = (2*k/(a*dx^2)+2*k/dy^2-rho*c*u(1,i)/dx)^(-1)*(k*((2/dx^2)*(T_d/(a*(a+1))+T(1,i+1)/(a+1))+2*T(2,i)/dy^2)-rho*c*u(1,i)*T_d/(a*dx));
               end
           end   
        end
    end
    
    %Margem superior
    for i = 2:size(X,2)-1
        if u(end,i) >= 0
            T(end,i) = (2*k/dx^2 + 2*k/dy^2 + rho*c*u(end,i))^(-1)*(k*((T(end,i+1)+T(end,i-1))/(dx^2) + 2*T(end-1,i)/(dy^2)) + rho*c*u(end,i)*T(end,i-1)/dx);
        else
            T(end,i) = (2*k/dx^2 + 2*k/dy^2 - rho*c*u(end,i))^(-1)*(k*((T(end,i+1)+T(end,i-1))/(dx^2) + 2*T(end-1,i)/(dy^2)) - rho*c*u(end,i)*T(end,i+1)/dx);
        end 
    end
               
    %Margem direita
    for j = 2:size(Y,1)-1
        T(j,end)=(2*k/dx^2 + 2*k/dy^2)^(-1)*(k*(2*T(j,end-1)/dx^2 + (T(j+1,end)+T(j-1,end))/dy^2));
    end
    
    %Canto superior direito
    T(end,end) = (2/dx^2 + 2/dy^2)^(-1) * (2*T(end, end-1)/dx^2 + 2*T(end-1,end)/dy^2);
    
    %Canto inferior direito
    T(1,end) = (2/dx^2 + 2/dy^2)^(-1) * (2*T(1, end-1)/dx^2 + 2*T(2,end)/dy^2);
end
