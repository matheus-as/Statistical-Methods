
%%Funcão de Ajuste Normal
function [Order1Sample,XTable] = NormalAdjustment(Sample)

X = Sample;
X(:,1)=1;
InitX = X;

ColNames = cell(1,1);
ColNames{1,1}=['x' num2str(0)];

for i=1:width(InitX)-1
  ColNames{height(ColNames)+1,1}=['x' num2str(i)];
end

XTable = array2table(X);
XTable.Properties.VariableNames = ColNames;
fprintf("<strong>Modelo de 1ª Ordem</strong>\n\n");
% disp(XTable);
Order1Sample=X;
Order1Sample(:,1) = Sample(:,1);
end