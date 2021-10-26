load fisheriris
x = meas(51:end,1:2);        % iris data, 2 classes and 2 features
y = (1:100)'>50;             % versicolor=0, virginica=1
b = glmfit(x,y,'binomial');  % logistic regression
p = glmval(b,x,'logit');     % get fitted probabilities for scores

figure
[X,Y] = perfcurve(species(51:end,:),p,'virginica');
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC for classification by logistic regression')

% Obtain errors on TPR by vertical averaging
[X,Y] = perfcurve(species(51:end,:),p,'virginica','nboot',10000,'xvals','all', 'BootType', 'normal');
errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1)); % plot errors