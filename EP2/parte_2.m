%--------------------------------------------------------------------------
% PMR3401 - SEGUNDO EXERCÍCIO PROGRAMA
% PARTE II - ITEM (a): ANÁLISE MODAL (FREQUÊNCIAS E MODOS DE VIBRAÇÃO)
%--------------------------------------------------------------------------
clc;
clear;
close all;

%% 1. PARÂMETROS E GEOMETRIA (CÓDIGO ANTERIOR)
% --- Parâmetros ---
m = 0.04; p = 0.40; t = 0.12; % NACA 4412 
N_pontos = 500; % Resolução de x: 500 valores. 

% --- Propriedades do Material e Seção ---
E = 210e9;      % Módulo de Young (Aço) [Pa] 
rho = 7800;     % Densidade (Aço) [kg/m^3] - Usado 7800 conforme código anterior 

% Pórtico (perfil)
b_portico = 0.05; % Lado da seção [m] 
A_portico = b_portico^2;
I_portico = (b_portico^4) / 12;

% Treliça (interno)
b_trelica = 0.03; % Lado da seção [m] 
A_trelica = b_trelica^2;
% I_trelica = (b_trelica^4) / 12; % Mesmo que não seja usado para treliça pura, não faz mal manter

% --- Geração de Nós ---
beta = linspace(0, pi, N_pontos+1);
x = (1 - cos(beta)) / 2; % Distribuição cossenoidal para x 

yt = 5 * t * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4); % Equação para yt(x) 

yc = zeros(size(x)); dyc_dx = zeros(size(x));
idx_p1 = find(x < p); 
idx_p2 = find(x >= p);

yc(idx_p1) = (m/p^2) * (2*p*x(idx_p1) - x(idx_p1).^2); % Equação para yc(x) se x < p 
dyc_dx(idx_p1) = (2*m/p^2) * (p - x(idx_p1)); % Equação para dyc/dx se x < p 

yc(idx_p2) = (m/(1-p)^2) * ((1-2*p) + 2*p*x(idx_p2) - x(idx_p2).^2); % Equação para yc(x) se x >= p 
dyc_dx(idx_p2) = (2*m/(1-p)^2) * (p - x(idx_p2)); % Equação para dyc/dx se x >= p 

theta = atan(dyc_dx); % Ângulo de inclinação da linha média 

xu = x - yt .* sin(theta); yu = yc + yt .* cos(theta); % Coordenadas do contorno superior 
xl = x + yt .* sin(theta); yl = yc - yt .* cos(theta); % Coordenadas do contorno inferior 

% Organização dos nós (primeiro os do contorno inferior, depois os do superior)
Nodes = zeros(2*(N_pontos+1), 3);
Nodes(1:N_pontos+1, :) = [(1:N_pontos+1)', fliplr(xl)', fliplr(yl)']; % Contorno inferior (do bordo de fuga para o de ataque)
Nodes(N_pontos+2:end, :) = [(N_pontos+2:2*(N_pontos+1))', xu', yu']; % Contorno superior (do bordo de ataque para o de fuga)

% --- Conectividade dos Elementos ---
% Elementos de Pórtico: Ligam nós adjacentes no mesmo contorno (superior ou inferior)
Elements_Portico = zeros(2*N_pontos, 3); % [ID, Nó1, Nó2]
for i = 1:N_pontos
    Elements_Portico(i, :) = [i, i, i+1]; % Elementos do contorno inferior
    Elements_Portico(N_pontos+i, :) = [N_pontos+i, N_pontos+1+i, N_pontos+1+i+1]; % Elementos do contorno superior
end

% Linhas guia para a estrutura interna 
x_guias = [0.05, 0.13, 0.17, 0.30, 0.40, 0.50, 0.60, 0.64, 0.75, 0.83];

% Encontrar os nós mais próximos das linhas guia nos contornos
nos_guia_inf = zeros(size(x_guias)); 
nos_guia_sup = zeros(size(x_guias));
for i = 1:length(x_guias)
    % Para o contorno inferior
    [~, idx_inf] = min(abs(Nodes(1:N_pontos+1, 2) - x_guias(i))); 
    nos_guia_inf(i) = Nodes(idx_inf, 1); % Guarda o ID do nó
    
    % Para o contorno superior (ajusta o índice para o bloco correto de nós)
    [~, idx_sup] = min(abs(Nodes(N_pontos+2:end, 2) - x_guias(i))); 
    nos_guia_sup(i) = Nodes(idx_sup + (N_pontos+1), 1); % Guarda o ID do nó
end

% Elementos de Treliça: Conectividade baseada na Figura 3 e nas linhas guia 
% Note que os índices em Elements_Trelica são os IDs dos nós, não os índices da matriz Nodes
Elements_Trelica = [
 1, nos_guia_inf(1), nos_guia_sup(2);
 2, nos_guia_inf(2), nos_guia_sup(2);
 3, nos_guia_sup(3), nos_guia_inf(3);
 4, nos_guia_sup(3), nos_guia_inf(4);
 5, nos_guia_inf(4), nos_guia_sup(5);
 6, nos_guia_sup(5), nos_guia_inf(6);
 7, nos_guia_inf(6), nos_guia_sup(7);
 8, nos_guia_inf(7), nos_guia_sup(7);
 9, nos_guia_sup(8), nos_guia_inf(8);
 10, nos_guia_inf(9), nos_guia_sup(8);
 11, nos_guia_sup(10), nos_guia_inf(9);
 12, nos_guia_sup(10), nos_guia_inf(10)
];

%% 2. MONTAGEM DAS MATRIZES GLOBAIS [K] E [M]
num_nodes = size(Nodes, 1); % Total de nós
num_dof = 3 * num_nodes; % 3 GDLs por nó (dx, dy, d_theta_z)

K_global = zeros(num_dof, num_dof); % Matriz de rigidez global
M_global = zeros(num_dof, num_dof); % Matriz de massa global

% --- Elementos de Pórtico --- 
for i = 1:size(Elements_Portico, 1)
    n1_idx = Elements_Portico(i, 2); % ID do primeiro nó
    n2_idx = Elements_Portico(i, 3); % ID do segundo nó
    
    n1_coords = Nodes(n1_idx, 2:3); % Coordenadas (x, y) do nó 1
    n2_coords = Nodes(n2_idx, 2:3); % Coordenadas (x, y) do nó 2
    
    dx = n2_coords(1) - n1_coords(1);
    dy = n2_coords(2) - n1_coords(2);
    L = sqrt(dx^2 + dy^2); % Comprimento do elemento
    
    c = dx/L; s = dy/L; % Cosseno e seno do ângulo de inclinação
    
    % Matrizes locais (K e M) de pórtico (elemento de viga de Euler-Bernoulli)
    EAL = E*A_portico/L; 
    EIL = E*I_portico/L; 
    EIL2 = E*I_portico/L^2; 
    EIL3 = E*I_portico/L^3;
    
    k_local = [
        EAL   0      0      -EAL  0      0;
        0     12*EIL3 6*EIL2 0     -12*EIL3 6*EIL2;
        0     6*EIL2  4*EIL  0     -6*EIL2  2*EIL;
        -EAL  0      0      EAL   0      0;
        0     -12*EIL3 -6*EIL2 0     12*EIL3 -6*EIL2;
        0     6*EIL2  2*EIL  0     -6*EIL2  4*EIL
    ];
    
    rhoAL = rho*A_portico*L/420; % Coeficiente para a matriz de massa consistente
    m_local = rhoAL * [
        140 0     0     70  0     0;
        0   156   22*L  0   54    -13*L;
        0   22*L  4*L^2 0   13*L  -3*L^2;
        70  0     0     140 0     0;
        0   54    13*L  0   156   -22*L;
        0   -13*L -3*L^2 0  -22*L 4*L^2
    ];
    
    % Matriz de Rotação (transforma coordenadas locais para globais)
    T = [c s 0 0 0 0; 
         -s c 0 0 0 0; 
         0 0 1 0 0 0; 
         0 0 0 c s 0; 
         0 0 0 -s c 0; 
         0 0 0 0 0 1];
    
    k_global_elem = T' * k_local * T;
    m_global_elem = T' * m_local * T;
    
    % Montagem na matriz global
    % Os GDLs para o nó n são (3n-2, 3n-1, 3n)
    dofs = [3*n1_idx-2, 3*n1_idx-1, 3*n1_idx, 3*n2_idx-2, 3*n2_idx-1, 3*n2_idx];
    K_global(dofs, dofs) = K_global(dofs, dofs) + k_global_elem;
    M_global(dofs, dofs) = M_global(dofs, dofs) + m_global_elem;
end

% --- Elementos de Treliça --- 
for i = 1:size(Elements_Trelica, 1)
    n1_idx = Elements_Trelica(i, 2); % ID do primeiro nó
    n2_idx = Elements_Trelica(i, 3); % ID do segundo nó
    
    n1_coords = Nodes(n1_idx, 2:3); % Coordenadas (x, y) do nó 1
    n2_coords = Nodes(n2_idx, 2:3); % Coordenadas (x, y) do nó 2
    
    dx = n2_coords(1) - n1_coords(1);
    dy = n2_coords(2) - n1_coords(2);
    L = sqrt(dx^2 + dy^2); % Comprimento do elemento
    
    c = dx/L; s = dy/L; % Cosseno e seno do ângulo de inclinação
    
    % Matrizes locais (K e M) de treliça (apenas GDLs de translação)
    EAL = E*A_trelica/L;
    k_local = EAL * [
        c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
        -c^2 -c*s c^2 c*s;
        -c*s -s^2 c*s s^2
    ];
    
    rhoAL = rho*A_trelica*L/6; % Coeficiente para a matriz de massa consistente de treliça
    m_local = rhoAL * [
        2 0 1 0;
        0 2 0 1;
        1 0 2 0;
        0 1 0 2
    ];
    
    % Montagem na matriz global (treliça só tem GDLs de translação: dx, dy)
    % Os GDLs para o nó n são (3n-2, 3n-1) para translação
    dofs = [3*n1_idx-2, 3*n1_idx-1, 3*n2_idx-2, 3*n2_idx-1];
    K_global(dofs, dofs) = K_global(dofs, dofs) + k_local;
    M_global(dofs, dofs) = M_global(dofs, dofs) + m_local;
end

%% 3. APLICAÇÃO DAS CONDIÇÕES DE CONTORNO
% Fixar a estrutura restringindo todos os graus de liberdade dos nós que
% pertencem à superfície superior do aerofólio e que estão compreendidos na
% região das longarinas. 
x_longarina1 = [0.13, 0.17]; % Intervalo da primeira longarina em x
x_longarina2 = [0.60, 0.64]; % Intervalo da segunda longarina em x

nos_restringidos_ids = []; % IDs dos nós a serem restringidos

% Percorre apenas os nós do contorno superior
nos_superiores_indices = N_pontos+2 : 2*(N_pontos+1); 
for i = 1:length(nos_superiores_indices)
    no_id = Nodes(nos_superiores_indices(i), 1); % ID do nó
    x_coord = Nodes(nos_superiores_indices(i), 2); % Coordenada x do nó
    
    % Verifica se o nó está nas regiões das longarinas
    if (x_coord >= x_longarina1(1) && x_coord <= x_longarina1(2)) || ...
       (x_coord >= x_longarina2(1) && x_coord <= x_longarina2(2))
        nos_restringidos_ids = [nos_restringidos_ids, no_id];
    end
end

% Converte IDs dos nós em GDLs a serem restringidos (dx, dy, d_theta_z para cada nó)
GDL_restringidos = [];
for i = 1:length(nos_restringidos_ids)
    no = nos_restringidos_ids(i);
    GDL_restringidos = [GDL_restringidos, 3*no-2, 3*no-1, 3*no];
end

% Define os GDLs livres (aqueles que não foram restringidos)
GDL_livres = setdiff(1:num_dof, GDL_restringidos);

% Reduz as matrizes globais K e M para considerar apenas os GDLs livres
K_red = K_global(GDL_livres, GDL_livres);
M_red = M_global(GDL_livres, GDL_livres);

%% 4. RESOLUÇÃO DO PROBLEMA DE AUTOVALOR
% Resolve o problema de autovalor generalizado: [K_red]{phi} = omega^2[M_red]{phi}
% Retorna os autovetores (Modos) e uma matriz diagonal de autovalores (Omega2)
[Modos, Omega2] = eig(K_red, M_red); 

% Extrai os autovalores (omega^2) da diagonal da matriz e os ordena em ordem crescente
omegas2 = diag(Omega2);
[omegas2_sorted, sort_idx] = sort(omegas2);
Modos_sorted = Modos(:, sort_idx); % Ordena os modos de vibração de acordo com as frequências

% Calcula as frequências naturais em rad/s e em Hz
omegas = sqrt(omegas2_sorted); % Frequências angulares (rad/s)
freqs_Hz = omegas / (2*pi); % Frequências em Hertz 

%% 5. APRESENTAÇÃO DOS RESULTADOS
fprintf('As 6 primeiras frequências naturais da estrutura são:\n'); 
for i = 1:6
    fprintf('  Modo %d: %.4f Hz\n', i, freqs_Hz(i));
end

% Plotagem dos modos de vibração 
% O código gerará uma figura para cada um dos 6 primeiros modos de vibração.
for i = 1:6
    figure('Name', sprintf('Modo de Vibração %d (%.2f Hz)', i, freqs_Hz(i)));
    hold on;
    
    % Desenha a estrutura não deformada em cinza (para referência)
    % Contorno do aerofólio (elementos de pórtico)
    for j = 1:size(Elements_Portico, 1)
        n1 = Elements_Portico(j,2); % ID do primeiro nó do elemento
        n2 = Elements_Portico(j,3); % ID do segundo nó do elemento
        plot([Nodes(n1,2), Nodes(n2,2)], [Nodes(n1,3), Nodes(n2,3)], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    % Estrutura interna (elementos de treliça)
    for j = 1:size(Elements_Trelica, 1)
        n1 = Elements_Trelica(j,2); % ID do primeiro nó do elemento
        n2 = Elements_Trelica(j,3); % ID do segundo nó do elemento
        plot([Nodes(n1,2), Nodes(n2,2)], [Nodes(n1,3), Nodes(n2,3)], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    
    % Monta o vetor de deslocamento completo (incluindo os GDLs restringidos como zero)
    U_completo = zeros(num_dof, 1);
    U_completo(GDL_livres) = Modos_sorted(:, i); % Preenche os GDLs livres com os valores do modo de vibração
    
    % Normaliza o modo de vibração para uma visualização consistente
    % Encontra o deslocamento máximo para normalização
    max_desloc = max(sqrt(U_completo(1:3:end).^2 + U_completo(2:3:end).^2));
    if max_desloc == 0 % Evita divisão por zero se o modo for todo zero (improvável para modos reais)
        escala = 0; 
    else
        escala = 0.1 / max_desloc; % Fator de escala para visualização (arbitrário, ajusta o tamanho da deformação)
    end

    % Calcula as novas coordenadas dos nós para a estrutura deformada
    % (deslocamentos em x e y são os 1o e 2o GDLs de cada nó)
    Nodes_deformados = Nodes(:, 2:3) + escala * [U_completo(1:3:end), U_completo(2:3:end)];
    
    % Desenha a estrutura deformada
    % Contorno do aerofólio (elementos de pórtico) em azul
    for j = 1:size(Elements_Portico, 1)
        n1 = Elements_Portico(j,2); 
        n2 = Elements_Portico(j,3);
        plot([Nodes_deformados(n1,1), Nodes_deformados(n2,1)], [Nodes_deformados(n1,2), Nodes_deformados(n2,2)], 'b-', 'LineWidth', 1.5);
    end
    % Estrutura interna (elementos de treliça) em vermelho
    for j = 1:size(Elements_Trelica, 1)
        n1 = Elements_Trelica(j,2); 
        n2 = Elements_Trelica(j,3);
        plot([Nodes_deformados(n1,1), Nodes_deformados(n2,1)], [Nodes_deformados(n1,2), Nodes_deformados(n2,2)], 'r-', 'LineWidth', 1.5);
    end
    
    title(sprintf('Modo de Vibração %d: %.2f Hz', i, freqs_Hz(i))); % Título com frequência 
    xlabel('Posição na Corda (x)'); % Rótulo do eixo X 
    ylabel('Posição Vertical (y)'); % Rótulo do eixo Y 
    axis equal; % Garante que os eixos tenham a mesma escala
    grid on; % Adiciona grade
    hold off;
end