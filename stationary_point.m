clear
clc

RegressionCoefficients= [1.3667;
       -0.0046;
        0.0179;
       -0.0625;
       -0.1050;
       -0.2059;
       -0.1647;
       -0.1434;
       -0.1622;
        0.1069;
        0.0419;
       -0.0731;
        0.0644;
       -0.0431;
       -0.1556;
       -0.0806;
        0.0144];


% RegressionCoefficients = [  79.9400;
%                             0.9951;
%                             0.5152;
%                            -1.3764;
%                            -1.0013;
%                             0.2500];
                        
[b,B,xs,ys] = stat_point(RegressionCoefficients);
printParameters(b,B,xs,ys);
CanonicalAnalyses(ys,B);

%% Calcula o Ponto Estacionário
function [b,B,xs,ys] = stat_point(RegressionCoefficients)
b0 = RegressionCoefficients(1);
RegressionCoefficients(1) = []; % Remove o primeiro coeficientes
dim = floor(sqrt(height(RegressionCoefficients))); % Dimensão da matriz quadrada
b = zeros(dim,1);   % Matriz dos Coeficientes de 1ª Ordem
B = zeros(dim,dim); % Matriz dos Coeficientes de 2ª Ordem
k = (numel(B)- height(B))/2; % Determina quantidade de elementos de iterações duplas
RegressionCoefficients(2*dim+k+1:length(RegressionCoefficients)) = []; %Remove elementos com mais de duas interações


% Obtem o vetor coluna com os coeficientes de 1st Ordem
for i=1:dim
    b(i,1) = RegressionCoefficients(i);
end

% Obtem a diagonal com os coeficientes de 2nd Ordem
for i=1:dim
    for j=1:dim 
        if i==j
            B(i,j)=RegressionCoefficients(dim+i);
        end
    end
end

v = RegressionCoefficients(2*dim+1:length(RegressionCoefficients));
v = (v')./2;

M=zeros(dim);       % Cria matriz M para preencher a região triangular inferior
triang_area = tril(~M,-1);   % Obtem a parte triangular inferior da matriz abaixo da diagonal principal 
B(triang_area) = v;          
B = B';             % Preenche a região triangular superior com os valores de v
B(triang_area)=v;   % Preenche a região triangular inferior com os valores de v

xs = -0.5*(B\b);    % calcula o ponto estacionário
ys = b0 +0.5*xs'*b; % calcula o y previsto

end

%% Imprime os parâmetros calculados
function printParameters(b,B,xs,ys)
fprintf("Coeficientes de 1ª ordem\n");
disp("b: ");
disp(b);
fprintf("Coeficientes de 2ª ordem (diagonal) e interação/2 (outra diagonal)\n");
disp("B: ");
disp(B);

fprintf("Ponto Estacionário\n");
fprintf("xs = ");
disp(xs');
fprintf("Valor de y previsto\n");
fprintf("ys = b0 +0.5*xs'*b = %f\n", ys);

end

%% Análise Canônica
function CanonicalAnalyses(ys,B)
eigen_values = flipud(eig(B));
fprintf("\nAs raízes desta equação são: ");
disp(eigen_values');
end