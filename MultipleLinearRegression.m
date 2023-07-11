%% Função Geral 
% Interface com usuário obtém os parametros de entrada para calcular a
% regressão multipla
function RegressionCoefficients = MultipleLinearRegression(Sample)
ENCODED_VARIABLES = 1;
NATURAL_VARIABLES = 2;
var = ENCODED_VARIABLES;

%% Nível de Significância
alpha=0.00;
%level_of_sig=1-alpha;
%critical_value=finv(level_of_sig,k,n-k-1);

while alpha <= 0.00
 prompt = {'Digite um valor para alpha:'};
 dlgtitle = 'Alfa';
 dims = [1 35];
 definput = {'0.01'};
 answer = inputdlg(prompt,dlgtitle,dims,definput);
 alpha = str2double(answer);
 fprintf("Valor de alpha: %g\n\n", alpha);
end

tic
%% Interface para verificar se será necessário realizar o ajuste quadrático da amostra
quest = 'Realizar ajuste quadrático sobre a amostra:';
answer = questdlg(quest, ....
	'Ajuste da Amostra', ....
	'Sim','Não','Não');
% Handle response
switch answer
    case 'Sim'
        [Order2Sample,XTable] = QuadraticAdjustment(Sample);
%         disp("Ajuste quadrático aplicado a amostra");
        fprintf("\n");
        previousSample = Sample;
        Sample = Order2Sample;
    case 'Não'
%         disp('Análise realizada sem o ajuste quadrático da amostra')
        previousSample = Sample;
        [Sample, XTable] = NormalAdjustment(Sample);
end
writetable(XTable,'D:\ArquivosMatlab\Estatística\Sample.txt'); % x is the name of your table, the csv will be saved where your code is saved
disp(XTable);

%% Calculo dos Parametros

% Sample = [ 0.87	1	1	1	1	1	1;
%            0.74	1	1	1	1	-1	1;
%            0.51	1	1	1	1	-1	1;
%            0.99	1	1	1	1	1	1;
%            0.67	1	1	1	1	1	-1;
%            0.72	1	1	1	1	-1	-1;
%            0.81	1	1	1	1	-1	-1;
%            1.01	1	1	1	1	1	-1;
%            1.33	1	1	1	1	1	-1;
%            0.7	1	1	1	1	-1	-1;
%            0.82	1	1	1	1	-1	-1;
%            0.78	1	1	1	1	1	-1;
%            0.36	1	1	1	1	1	1;
%            0.23	1	1	1	1	-1	1;
%            0.21	1	1	1	1	-1	1;
%            0.44	1	1	1	1	1	1;
%            0.56	4	0	0	0	0	0;
%            0.49	4	0	0	0	0	0;
%            0.57	0	4	0	0	0	0;
%            0.81	0	4	0	0	0	0;
%            0.9	0	0	4	0	0	0;
%            0.65	0	0	4	0	0	0;
%            0.91	0	0	0	4	0	0;
%            0.49	0	0	0	4	0	0;
%            1.43	0	0	0	0	0	0;
%            1.17	0	0	0	0	0	0;
%            1.50	0	0	0	0	0	0];
       
n = getObservations(Sample);
k = getRegressionDoF(Sample);
DoFResidue = getDoFResidue(n,k);

FTable = getFTable(Sample,n,k);
RegressionTable = getRegresionTable(Sample,XTable,n,k);

fprintf('<strong>Tabela ANOVA</strong>\n\n');
disp(FTable);
fprintf("Para alpha = %g, verifica-se que %s", alpha,getANOVAResult(FTable,alpha));

fprintf('\n\n<strong>Tabela dos Coeficientes</strong>\n\n');
disp(RegressionTable);

R2 = fR2(FTable);
R2Adj = fR2Adj(FTable);

ResiduesTable = getResiduesTable(Sample,RegressionTable);

%% Interface para verificar se deve imprimir a Análise dos Resíduos
quest = 'Deseja imprimir os gráficos da análise dos resíduos:';
answer = questdlg(quest, ....
	'Gráficos dos Resíduos', ....
	'Sim','Não','Não');
% Handle response
switch answer
    case 'Sim'
        quest = 'Os dados foram coletados de maneira ordenada?';
        answer = questdlg(quest, ....
        'ordem de coleta dos dados', ....
        'Sim','Não','Não');
        printResidueGraphs(previousSample,var,ResiduesTable, answer);
    case 'Não'
end

printCoeficientsResult(Sample,RegressionTable, alpha,n,k);

RegressionCoefficients = RegressionTable.Beta;

RegressionResults(Sample, FTable, R2, R2Adj, n,k);
toc
end

%% Calculo de n
% Obter o número de observações
function n = getObservations(Sample)
 X = Sample;
 X(:,1) = 1;
 
 n = height(X);      %número de amostras
end

%% Calculo de k [DoF da Regressao]
% Obtem o numero de graus de liberdade da amostra
function k = getRegressionDoF(Sample)
 X = Sample;
 X(:,1) = 1;

 k = width(X)-1;     %número de variáveis
end

%% Calculo de DoF Residuos
% Obtem o numero de graus de liberdade dos resíduos
function DoFResidue = getDoFResidue(n,k)
 DoFResidue = n - k - 1;
end

%% Calcula a ANOVA 
% Obtem uma tabela com a análise de variância
function FTable = getFTable(Sample,n,k)
 X = Sample;
 X(:,1) = 1; 
 y = Sample(:,1);
 Hat = (inv(X'*X))*X';
 beta = Hat*y;
 
% Calcula o Somatório de Y
 y_total = 0;
 for i=1:length(y)
  y_total=y_total+y(i);
 end

 SSregressao = beta'*X'*y - ((y_total)^2)/length(y); % Soma dos Quadrados da Regressão
 SStotal = y'*y - ((y_total)^2)/length(y);           % Soma dos Quadrados Totais
 SSres = y'*y - beta'*X'*y;                          % Soma dos Resíduos

 F = ((SSregressao/k)/(SSres/(n-k-1)));  % valor F
 p_value = 1 - fcdf(F,k,n-k-1);          % p-valor
 
 % Tabela da ANOVA
 Data = [  k         SSregressao     SSregressao/k       F   p_value
           n-k-1      SSres           SSres/(n-k-1)      0   0
           n-1       SStotal         0                   0   0];
 VarNames = {'Graus de Liberdade', 'Soma de Quadrados', 'Quadrado Medio', 'F', 'pValor'};
 SourceOfVariation = {'Regressao','Residuo','Total'};
 FTable = table(Data(:,1),Data(:,2),Data(:,3),Data(:,4), Data(:,5), 'VariableNames',VarNames,'RowNames',SourceOfVariation);
end

%% Verifica se o modelo é significativo por meio da ANOVA
% Interpreta os resultados a partir da Tabela F que foi calculada pela
% função getFTable
function ANOVAResult = getANOVAResult(FTable,alpha)
 FTable.Properties.VariableNames = {'DoF', 'SumOfSquares', 'MeanSquare','F', 'pValor'};
 p_value = FTable.pValor(1);
 F = FTable.F(1);
 
 if p_value<alpha
  ANOVAResult = "as variáveis consideradas são estatiscamente significativas para o modelo," + newline + "pois P-Valor("+F+")= " + p_value+ "<" +alpha;
 elseif  p_value>alpha
  ANOVAResult = "as variáveis consideradas não são estatiscamente significativas para o modelo," + newline+" pois P-Valor("+F+")= " + p_value+ ">" +alpha;
 end
end

%% Calculo a Tabela dos Coeficientes de Regressao
% Obtem uma tabela com os parâmetros da regressão
function RegressionTable = getRegresionTable(Sample,XTable,n,k)
 X = Sample;
 X(:,1) = 1;
 y = Sample(:,1);
 
 XRowNames = XTable.Properties.VariableNames; %Obtem o nome das variáveis

 tb=[width(X),1];         % Vetor t-Student
 erropadrao=[width(X),1];
 C = inv(X'*X);           % Obtem a matriz C
 Hat = (inv(X'*X))*X';
 beta = Hat*y;
 
 SSres = y'*y - beta'*X'*y; % Calculo da Soma dos Resíduos
 MSRes = SSres/(n-k-1);   % Quadrado Medio Residual

 for i=1:length(beta)
  tb(i)= beta(i)/sqrt(MSRes*C(i,i));
  erropadrao(i)=sqrt(MSRes*C(i,i));
 end
 
tb=tb';
erropadrao = erropadrao';

abs_tb = abs(tb);        % Absolute values t-student values
cp = tcdf(abs_tb,n-k-1); % Computes Cumulative Probability
pValue = 2*(1 - cp);     % Calcula os pvalores
RowNames = cell(width(X)-1,1);
% Preenche coluna de Variáveis do Teste dos Coeficientes da Regressão
 for i=1:(width(X))
  if i==1
   RowNames{i,1}='(Intercept)';
  else
   RowNames{i,1}=char(XRowNames(i));
  end
 end

VarNames = {'Beta', 'Erro Padrao', 'Stat t', 'valor-P'};
RegressionTable = table(beta,erropadrao,tb,pValue, 'VariableNames',VarNames,'RowNames',RowNames);
end

%% Coeficient Analyses
% Interpreta os resultados da regressão da Tabela dos Coeficientes 
function printCoeficientsResult(Sample,RegressionTable,alpha,n,k)
X = Sample;

RegressionTable.Properties.VariableNames = {'Beta', 'ErroPadrao', 't', 'valorP'};
XRowNames = RegressionTable.Properties.RowNames;
validBetas = [];

% disp(char(XRowNames(1)));

fprintf(" Assim a partir da tabela dos Coeficientes da Regressão, para alpha = %g, verifica-se que:\n\n", alpha);

% Verifica Variáveis significativas para o Modelo
for i=1:length(RegressionTable.Beta)
    if RegressionTable.valorP(i) < alpha
     fprintf("  B[%d] = %.4f é estatiscamente significativa para o modelo, pois P-Valor(%.4g,%d) = %.4f < %.4f.\n",i-1, RegressionTable.Beta(i),RegressionTable.t(i),n-k-1,RegressionTable.valorP(i),alpha);
     validBetas(width(validBetas)+1) = RegressionTable.Beta(i);
    end
    if i==1
       XRowNames(i) = {''};
    end
    if RegressionTable.valorP(i)>alpha
     fprintf("  B[%d] = %.4f não é estatiscamente significativa para o modelo, pois P-Valor(%.4g,%d) = %.4f > %.4f.\n",i-1, RegressionTable.Beta(i),RegressionTable.t(i),n-k-1,RegressionTable.valorP(i),alpha);    
     XRowNames(i) = {''};
    end
end
 
 XRowNames = XRowNames(~cellfun('isempty',XRowNames)); % identify and remove empty cells
%  disp(validBetas);
 
  if isempty(XRowNames) %verifica se existem variavies x
    fprintf("\n Dessa forma, verifica-se que nenhuma das variáveis representa o modelo.\n");
  else 
   for i=1:length(XRowNames)
      if i==1
        fprintf("\nDessa forma, verifica-se que o modelo deve conter apenas: ");
        fprintf("%s, ", char(XRowNames(i)));
      elseif i>1 && i<=(length(XRowNames)-1)
        fprintf("%s, ", char(XRowNames(i)));
      elseif i==length(XRowNames)
        fprintf("e %s. \n", char(XRowNames(i)));
      end
   end
   fprintf("Com isso, tem-se modelo sugerido: ");
   for i=1:length(validBetas)
        if i == 1
            fprintf("y = %.4f", (validBetas(i)));
        end
        if i>1 
            if validBetas(i) > 0
                fprintf("+%.4f%s", validBetas(i),char(XRowNames(i-1)));
            elseif validBetas(i)<0
                fprintf("%.4f%s", validBetas(i),char(XRowNames(i-1)));
            end
        end
        if i==(length(validBetas))
                fprintf(".\n\n");
        end
   end
  end
end

%% Obtem a Tabela dos Residuos
% Calcula a Tabela dos Residuos a partir da Tabela de Regressão
function ResiduesTable = getResiduesTable(Sample,RegressionTable)
 X = Sample;
 X(:,1) = 1;      % obtem os valores da variavel X
 y = Sample(:,1); % obtem os valores de y

 RegressionTable.Properties.VariableNames = {'Beta', 'ErroPadrao', 't', 'valorP'};
 
 yest = X*(RegressionTable.Beta);
 e = y-yest;
 e_ord = sort(e);

 p=(length(y):1);   % vetor de probabilidade           
 pos=(length(y):1); % vetor de posição

 for i=1:length(y)   
  p(i)=((i-0.5)*100)/length(y);
  pos(i) = i;
 end
 
 p=p';
 pos=pos';

 DataResidue = [y yest e pos e_ord p];
 VarNames = {'y', 'Valor Ajustado', 'Resíduo', 'i', 'Resíduo Ordenado', 'p(%)'};
 ResiduesTable = table(DataResidue(:,1),DataResidue(:,2),DataResidue(:,3),DataResidue(:,4), DataResidue(:,5), DataResidue(:,6), 'VariableNames',VarNames);
end

%% Imprime Gráficos dos Resíduos
function printResidueGraphs(Sample,var,ResiduesTable, answer)

 X = Sample;
 X(:,1) = 1;
 % Ajusta o nome das variáveis
 VarNames = {'y', 'yest', 'residuo', 'i', 'e_ord', 'p'};
 ResiduesTable.Properties.VariableNames = VarNames;
 
 ENCODED_VARIABLES = 1;
%  NORMAL_VARIABLE = 2;
 
 figure;
 hold on;
 title('NPP dos Resíduos');
 scatter(norminv((ResiduesTable.p)/100),ResiduesTable.e_ord,'filled');
 title('NPP dos Resíduos');
 xlabel('Z normal');
 ylabel('Resíduos ordenados');
 hold off;
 
 figure;
 hold on;
 title('NPP dos Resíduos');
 scatter(ResiduesTable.p,ResiduesTable.e_ord,'filled');
 title('NPP dos Resíduos');
 xlabel('Probabilidade');
 ylabel('Resíduos ordenados');
 hold off;

 % A partir do Gráfico de Distribuição Normal dos Resíduos verifica-se que a distribuição
 %plot Residuals vs Fitted Values 
 figure
 title('Análise de Normalidade');
 scatter(ResiduesTable.yest,ResiduesTable.residuo,'filled');
 title('Resíduos vs Valores Ajustados');
 xlabel('Valores ajustados');
 ylabel('Resíduos');
 
 if var == ENCODED_VARIABLES
   for i=1:width(X)-1
     figure
     title('Analise da Variancia dos Residuos');
     scatter(X(:,i+1),ResiduesTable.residuo,'filled');
     
     str_x = ['X' num2str(i)];
     GraphTitle = ['Resíduos vs ' str_x];
     title(GraphTitle);
     xlabel(str_x);
     ylabel('Resíduos');
   end
 end
 index = (X(:,1));
 switch (answer)
     case 'Sim'
         for i=1:height(X)
             index(i) = i;
         end
         figure
         title('Analise da Independência dos Dados');
         scatter(index,ResiduesTable.residuo,'filled');
         GraphTitle = ['Resíduos vs Ordem de Coleta'];
         title(GraphTitle);
         xlabel('Ordem');
         ylabel('Resíduos');
     case 'Não'
 end

end

%% Calculo de R Quadrado
function R2 = fR2(FTable)
 FTable.Properties.VariableNames = {'DoF', 'SumOfSquares', 'MeanSquare','F', 'pValor'};
 
 SSres = FTable.SumOfSquares(2);
 SStotal = FTable.SumOfSquares(3);
 
 R2 = 1 - SSres/SStotal;
end

%% Calculo de R Quadrado Adjustado
function R2Adj = fR2Adj(FTable)
 FTable.Properties.VariableNames = {'DoF', 'SumOfSquares', 'MeanSquare','F', 'pValor'};

 SSres = FTable.SumOfSquares(2);    % Soma dos Quadrados dos Residuos
 SStotal = FTable.SumOfSquares(3);  % Soma dos Quadrados Totais
 DoF_Residue = FTable.DoF(2);       % Graus de liberdade dos Resíduos
 DoF_Total= FTable.DoF(3);          % Graus de Liberdade Totais
 
 R2Adj = 1 - (SSres/DoF_Residue)/(SStotal/DoF_Total);
end

%% Imprime os resultados da regressão
function RegressionResults(Sample,FTable,R2,R2Adj, n,k)
 FTable.Properties.VariableNames = {'DoF', 'SumOfSquares', 'MeanSquare','F', 'pValor'};
 y = Sample(:,1);
 SSres = FTable.SumOfSquares(2);
 F = FTable.F(1);
 DoF_Regression = FTable.DoF(1);
 DoF_Residue = FTable.DoF(2);
 pValor = FTable.pValor(1);
 
 fprintf("<strong>Resultados da Regressão</strong>\n\n");
 fprintf("Método: Mínimos Quadrados");
 fprintf("\nF-statistic: %f", F);
 fprintf("\nProb (F-statistic): %g\n", pValor);
 fprintf("Número de observações = %d\n", length(y));
 fprintf("Graus de Liberdade dos Resíduos: %d\n",DoF_Residue);
 fprintf("Graus de Liberadade da Regressão: %d\n",DoF_Regression);
 fprintf("Erro Padrão = Raiz(MSRes)= %.5g\n", sqrt(SSres/(n-k-1)));
 fprintf("R2: %g", R2);
 fprintf("\nR2 Ajustado = %.4g, indicando que %.4g por cento da variabilidade dos dados é explicada pelo modelo.\n\n", R2Adj, R2Adj*100);
end
