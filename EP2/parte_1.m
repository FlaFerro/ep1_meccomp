clc, clear  

%% Plot do perfil

% Parâmetros NACA 4412
m = 0.04;   % Camber máximo (4%)
p = 0.4;    % Posição do camber máximo (40% da corda)
t = 0.12;   % Espessura máxima (12%)

% Gerar pontos ao longo da corda
beta = linspace(0, pi, 500);
x = (1 - cos(beta)) / 2;

% Distribuição de espessura
y_t = 5 * t * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);

% Linha média (camber) e sua derivada
mask = x < p;
y_c = zeros(size(x));
dy_cdx = zeros(size(x));

y_c(mask) = (m / p^2) * (2*p*x(mask) - x(mask).^2);
dy_cdx(mask) = (2*m / p^2) * (p - x(mask));

y_c(~mask) = (m / (1 - p)^2) * ((1 - 2*p) + 2*p*x(~mask) - x(~mask).^2);
dy_cdx(~mask) = (2*m / (1 - p)^2) * (p - x(~mask));

% Ângulo da linha média
theta = atan(dy_cdx);

% Coordenadas dos contornos
x_u = x - y_t .* sin(theta);
y_u = y_c + y_t .* cos(theta);
x_l = x + y_t .* sin(theta);
y_l = y_c - y_t .* cos(theta);

% Plotagem
figure;
hold on;
plot(x_u, y_u, 'r', 'LineWidth', 1.5, 'DisplayName', 'Contorno Superior');
plot(x_l, y_l, 'b', 'LineWidth', 1.5, 'DisplayName', 'Contorno Inferior');
plot(x, y_c, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Linha Média');
axis equal;
xlabel('Coordenada x normalizada');
ylabel('Coordenada y normalizada');
title('Perfil NACA 4412');
legend('Location', 'Northwest');
grid on;
hold off;

% Plot Cp

% Parâmetros adicionais
alpha = deg2rad(8);     % Ângulo de ataque em radianos
V = 300 * (1000/3600);  % Velocidade em m/s (300 km/h → 83.333 m/s)
rho = 1.29;             % Densidade do ar (kg/m³)
p_atm = 1e5;            % Pressão atmosférica (Pa)
q = 0.5 * rho * V^2;    % Pressão dinâmica

% Definir máscaras (mesmo código)
mask_2 = (x > 0.03) & (x < 0.97);

% Contorno Superior (Cp negativo)
Cp_upper = zeros(size(x));
Cp_upper(mask_2) = -2 * (alpha ./ sqrt(x(mask_2)) + dy_cdx(mask_2));

% Contorno Inferior (Cp positivo)
Cp_lower = zeros(size(x));
Cp_lower(mask_2) = 2 * (alpha ./ sqrt(x(mask_2)) + dy_cdx(mask_2));

% Plotagem 
figure;
hold on;
plot(x, Cp_upper, 'r', 'LineWidth', 1.5);
plot(x, Cp_lower, 'b', 'LineWidth', 1.5);
xlabel('x'); ylabel('C_p'); 
title('Coeficiente de Pressão - NACA 4412');
legend('Superior', 'Inferior');
grid on;

%% Pressao absoluta

p_lower = Cp_lower*q+p_atm;
p_upper = Cp_upper*q+p_atm;

figure;
hold on;
plot(x, p_upper, 'r', 'LineWidth', 1.5);
plot(x, p_lower, 'b', 'LineWidth', 1.5);
xlabel('x'); ylabel('Pressão (Pa)'); 
title('Pressão absoluta');
legend('Superior', 'Inferior');
grid on;

%% Forcas Fx e Fy

f_eq_lower = zeros(length(p_lower),2);
f_eq_upper = zeros(length(p_upper),2);

for i=1:length(p_lower)-1
    
    p_eq_lower = (p_lower(i)+p_lower(i+1))/2-p_atm;
    p_eq_upper = (p_upper(i)+p_upper(i+1))/2-p_atm;
    
    f_eq_lower(i,1)=f_eq_lower(i,1) - (p_eq_lower/2)*(y_l(i+1)-y_l(i));
    f_eq_lower(i+1,1)=f_eq_lower(i+1,1) - (p_eq_lower/2)*(y_l(i+1)-y_l(i));
    f_eq_lower(i,2)=f_eq_lower(i,2) + (p_eq_lower/2)*(x_l(i+1)-x_l(i));
    f_eq_lower(i+1,2)=f_eq_lower(i+1,2) + (p_eq_lower/2)*(x_l(i+1)-x_l(i));
    
    f_eq_upper(i,1)=f_eq_upper(i,1) + (p_eq_upper/2)*(y_u(i+1)-y_u(i));
    f_eq_upper(i+1,1)=f_eq_upper(i+1,1) + (p_eq_upper/2)*(y_u(i+1)-y_u(i));
    f_eq_upper(i,2)=f_eq_upper(i,2) - (p_eq_upper/2)*(x_u(i+1)-x_u(i));
    f_eq_upper(i+1,2)=f_eq_upper(i+1,2) - (p_eq_upper/2)*(x_u(i+1)-x_u(i));
    
end

% Parâmetros para visualização
N = 10;  % Usar 1 a cada N vetores
max_force = max([max(vecnorm(f_eq_upper, 2, 2)), max(vecnorm(f_eq_lower, 2, 2))]);
scale_factor = 0.1 / max_force;  % Menor escala para não exagerar no tamanho

% Indices espaçados para visualização
idx_u = 1:N:length(x_u);
idx_l = 1:N:length(x_l);

% Plot
figure;
hold on;

plot(x_u, y_u, 'k', 'LineWidth', 1.5);
plot(x_l, y_l, 'k', 'LineWidth', 1.5);

quiver(x_u(idx_u)', y_u(idx_u)', ...
       f_eq_upper(idx_u,1)*scale_factor, f_eq_upper(idx_u,2)*scale_factor, ...
       0, 'r', 'LineWidth', 1, 'MaxHeadSize', 0.3);

quiver(x_l(idx_l)', y_l(idx_l)', ...
       f_eq_lower(idx_l,1)*scale_factor, f_eq_lower(idx_l,2)*scale_factor, ...
       0, 'b', 'LineWidth', 1, 'MaxHeadSize', 0.3);

axis equal;
xlabel('Coordenada x normalizada');
ylabel('Coordenada y normalizada');
title('Vetores de força ao longo do perfil');
grid on;
hold off;
    
%% Forca de arrasto e sustentacao

T = [cos(alpha), sin(alpha);
    -sin(alpha), cos(alpha)];

f = [sum(f_eq_lower(:,1))+sum(f_eq_upper(:,1));
    sum(f_eq_lower(:,2))+sum(f_eq_upper(:,2))];

F = T*f;
    




