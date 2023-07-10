% Sample = [ 39.3  -1    -1;
%            40.0  -1     1;     
%            40.9   1    -1; 
%            41.5   1     1;  
%            40.3   0     0;  
%            40.5   0     0;  
%            40.7   0     0;  
%            40.2   0     0;  
%            40.6   0     0];

% Sample = [  66    -1    -1    -1    1;
%             70    -1    -1    1     1;
%             78    -1    1     -1    1;
%             60    -1    1     1     1;
%             80    1     -1    -1    -1;
%             70    1     -1    1     -1;
%             100   1     1     -1    -1;
%             75    1     1     1     -1;
%             100   -1.682      0     0     0;
%             80    1.682 0     0     0;
%             68    0     -1.682      0     0;
%             63    0     1.682 0     0;
%             65    0     0     -1.682      0;
%             82    0     0     1.682 0;
%             113   0     0     0     -1.682;
%             100   0     0     0     1.682;
%             118   0     0     0     0;
%             88    0     0     0     0;
%             100   0     0     0     0;
%             85    0     0     0     0];

% Sample = [66    -1    -1    -1;
%           70    -1    -1    1;
%           78    -1    1     -1;
%           60    -1    1     1;
%           80    1     -1    -1;
%           70    1     -1    1;
%           100   1     1     -1;
%           75    1     1     1;
%           100   -1.682      0     0;
%           80    1.682 0     0;
%           68    0     -1.682      0;
%           63    0     1.682 0;
%           65    0     0     -1.682;
%           82    0     0     1.682;
%           113   0     0     0;
%           100   0     0     0;
%           118   0     0     0;
%           88    0     0     0;
%           100   0     0     0;
%           85    0     0     0];


% Sample = [76.5	-1		 -1;
% 		  77	-1		 1;
% 		  78	 1		 -1;
% 		  79.5	 1		 1;
% 		  79.9	 0		 0;
% 		  80.3	 0		 0;
% 		  80	 0		 0;
% 		  79.7	 0		 0;
% 		  79.8	 0		 0;
% 		  78.4	 1.414	 0;
% 		  75.6	-1.414	 0;
% 		  78.5	 0	 1.414;
% 		  77	 0	-1.414];

% Sample =[66	-1	-1	-1;
% 		 70	-1	-1	1;
% 		 78	-1	1	-1;
% 		 60	-1	1	1;
% 		 80	1	-1	-1;
% 		 70	1	-1	1;
% 		 100	1	1	-1;
% 		 75     1	1	1;
% 		 100	-1.682	0	0;
% 		 80     1.682	0	0;
% 		 68     0	-1.682	0;
% 		 63     0	1.682	0;
% 		 65     0	0	-1.682;
% 		 82     0	0	1.682;
% 		 113	0	0	0;
% 		 100	0	0	0;
% 		 118	0	0	0;
% 		 88     0	0	0;
% 		 100	0	0	0;
% 		 85     0	0	0];

function [Order2Sample,XTable] = QuadraticAdjustment(Sample)

X = Sample;
X(:,1)=1;
InitX = X;

% NumOfArrays = width(X);
NumOfIterations = width(InitX)-2; % número de iterações duplas entre variáveis diferentes
ColNames = cell(1,1);
ColNames{1,1}=['x' num2str(0)];


% disp(NumOfArrays);
% disp(NumOfIterations); 

 for i=1:width(InitX)-1
    ColNames{height(ColNames)+1,1}=['x' num2str(i)];
%     disp(X);
 end

 for i=1:nchoosek(width(InitX)-1,1)
%     fprintf("iteração[%d]:x[%d]*x[%d] \n",i, i,i);
    X(:,width(X)+1)= InitX(:,i+1).^2;
    ColNames{height(ColNames)+1,1}=['x' num2str(i) '^2'];
%     disp(X);
  end
  
%   disp(NumOfArrays);
  
 pos = 1;
 for i=1:NumOfIterations
    for j=i:NumOfIterations
      if i<=nchoosek(width(InitX)-1,2)
%         fprintf("iteração[%d]:x[%d]*x[%d] \n",i,i,j+1);
        X(:,width(X)+pos) = InitX(:,i+1).*InitX(:,j+2);
        ColNames{height(ColNames)+1,1}=['x' num2str(i) '*x' num2str(j+1)];
%         disp(X);
      else
        pos = pos+1;
        break;
      end
    end
  end
  
% disp(ColNames);
if width(InitX)-1 >= 3
 i =1;
 NumOfGroups = nchoosek(width(InitX)-1,3);
 while i <= NumOfGroups
   LastPosition = i+3;
     if LastPosition <= width(InitX)
%       fprintf("iteração[%d]:x[%d]*x[%d]*x[%d] \n",i, i,i+1,i+2);
      X(:,width(X)+1) = InitX(:,i+2).*InitX(:,i+2).*InitX(:,i+3);
%       disp(X); 
      ColNames{height(ColNames)+1,1}=['x' num2str(i) '*x' num2str(i+1) '*x' num2str(i+2)];
     end
   i = i+1;
 end
end

XTable = array2table(X);
XTable.Properties.VariableNames = ColNames;
fprintf("<strong>Modelo de 2ª Ordem</strong>\n\n");
% disp(XTable);
Order2Sample=X;
Order2Sample(:,1) = Sample(:,1);

end