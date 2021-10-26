% 10/21/2021 Makoto. Cohen's d added for revision.
% 03/01/2021 Makoto. Used again.
% 02/08/2021 Makoto. Finished. Shaun's Excel data are screwed up. Fixed.
% 02/05/2021 Makoto. Created.

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract age-matched 40 patients and 20 controls. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NUM1,TXT1,RAW1] = xlsread('/data/mobi/Hiroki/code/QHYPS_020221_1647_exportMM.xlsx');

case40Idx    = find(contains(TXT1(1,:),'CASE40'));
control20Idx = find(contains(TXT1(1,:),'CONTROL20'));
preWake1Idx  = find(contains(TXT1(1,:),'CodePreWake1'));
preWake2Idx  = find(contains(TXT1(1,:),'CodePreWake2'));
preSleep1Idx = find(contains(TXT1(1,:),'CodePreSleep1'));
preSleep2Idx = find(contains(TXT1(1,:),'CodePreSleep2'));

case40ExcelSubjIdx      = find(NUM1(:,case40Idx-1)); % -1 is to adjust the difference between NUM and RAW.
control20ExcelSubjIdx   = find(NUM1(:,control20Idx-1));

preWake1ExcelCaseNames  = TXT1(case40ExcelSubjIdx+1, preWake1Idx);
preWake2ExcelCaseNames  = TXT1(case40ExcelSubjIdx+1, preWake2Idx);
preSleep1ExcelCaseNames = TXT1(case40ExcelSubjIdx+1, preSleep1Idx);
preSleep2ExcelCaseNames = TXT1(case40ExcelSubjIdx+1, preSleep2Idx);

preWake1ExcelControlNames  = TXT1(control20ExcelSubjIdx+1, preWake1Idx);
preWake2ExcelControlNames  = TXT1(control20ExcelSubjIdx+1, preWake2Idx);
preSleep1ExcelControlNames = TXT1(control20ExcelSubjIdx+1, preSleep1Idx);
preSleep2ExcelControlNames = TXT1(control20ExcelSubjIdx+1, preSleep2Idx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load precomputed MIs for Cleaning 0. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NUM2,TXT2,RAW2] = xlsread('/data/mobi/Hiroki/p4100_groupLevelStructure/rawData.xlsx');

    % Rename the files (Done only once)
    rawSubjNameList = TXT2(1,2:end)';
    renamingIdx = find(contains(rawSubjNameList, 'DQB93')); rawSubjNameList{renamingIdx} = 'DQB93a';
    renamingIdx = find(contains(rawSubjNameList, 'efv94')); rawSubjNameList{renamingIdx} = 'EFV94';
    renamingIdx = find(contains(rawSubjNameList, 'QGC43')); rawSubjNameList{renamingIdx} = 'QGC43a';
    renamingIdx = find(contains(rawSubjNameList, 'twg47')); rawSubjNameList{renamingIdx} = 'TWG47';
    renamingIdx = find(contains(rawSubjNameList, 'XBI94')); rawSubjNameList{renamingIdx} = 'XBI94a';
    renamingIdx = find(contains(rawSubjNameList, 'JDG25')); rawSubjNameList{renamingIdx} = 'JDG25a';
    renamingIdx = find(contains(rawSubjNameList, 'STV75')); rawSubjNameList{renamingIdx} = 'STV75a';
    renamingIdx = find(contains(rawSubjNameList, 'EGD95')); rawSubjNameList{renamingIdx} = 'EGD95a';
    renamingIdx = find(contains(rawSubjNameList, 'USP38')); rawSubjNameList{renamingIdx} = 'USP38a';
    renamingIdx = find(contains(rawSubjNameList, 'DEV35')); rawSubjNameList{renamingIdx} = 'DEV35a';
    renamingIdx = find(contains(rawSubjNameList, 'TDS79')); rawSubjNameList{renamingIdx} = 'TDS79a';

    [~,caseIdx1] = intersect(rawSubjNameList, preWake1ExcelCaseNames);
    [~,caseIdx2] = intersect(rawSubjNameList, preWake2ExcelCaseNames);
    [~,caseIdx3] = intersect(rawSubjNameList, preSleep1ExcelCaseNames);
    [~,caseIdx4] = intersect(rawSubjNameList, preSleep2ExcelCaseNames);
    [~,controlIdx1] = intersect(rawSubjNameList, preWake1ExcelControlNames);
    [~,controlIdx2] = intersect(rawSubjNameList, preWake2ExcelControlNames);
    [~,controlIdx3] = intersect(rawSubjNameList, preSleep1ExcelControlNames);
    [~,controlIdx4] = intersect(rawSubjNameList, preSleep2ExcelControlNames);

caseWake1Clean0     = NUM2(end,caseIdx1);
caseWake2Clean0     = NUM2(end,caseIdx2);
caseSleep1Clean0    = NUM2(end,caseIdx3);
caseSleep2Clean0    = NUM2(end,caseIdx4);
controlWake1Clean0  = NUM2(end,controlIdx1);
controlWake2Clean0  = NUM2(end,controlIdx2);
controlSleep1Clean0 = NUM2(end,controlIdx3);
controlSleep2Clean0 = NUM2(end,controlIdx4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load precomputed MIs for Cleaning 1. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUM3 = xlsread('/data/mobi/Hiroki/p4100_groupLevelStructure/clean1Data.xlsx');
caseWake1Clean1     = NUM3(end,caseIdx1);
caseWake2Clean1     = NUM3(end,caseIdx2);
caseSleep1Clean1    = NUM3(end,caseIdx3);
caseSleep2Clean1    = NUM3(end,caseIdx4);
controlWake1Clean1  = NUM3(end,controlIdx1);
controlWake2Clean1  = NUM3(end,controlIdx2);
controlSleep1Clean1 = NUM3(end,controlIdx3);
controlSleep2Clean1 = NUM3(end,controlIdx4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load precomputed MIs for Cleaning 2. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUM4 = xlsread('/data/mobi/Hiroki/p4100_groupLevelStructure/clean2Data.xlsx');
caseWake1Clean2     = NUM4(end,caseIdx1);
caseWake2Clean2     = NUM4(end,caseIdx2);
caseSleep1Clean2    = NUM4(end,caseIdx3);
caseSleep2Clean2    = NUM4(end,caseIdx4);
controlWake1Clean2  = NUM4(end,controlIdx1);
controlWake2Clean2  = NUM4(end,controlIdx2);
controlSleep1Clean2 = NUM4(end,controlIdx3);
controlSleep2Clean2 = NUM4(end,controlIdx4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Average Time 1 and Time 2. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace the outlier with the mean value.
caseWake2Clean0(25) = mean(caseWake2Clean0([1:24 26:end]));

caseWakeClean0 = (caseWake1Clean0 + caseWake2Clean0)/2;
caseWakeClean1 = (caseWake1Clean1 + caseWake2Clean1)/2;
caseWakeClean2 = (caseWake1Clean2 + caseWake2Clean2)/2;
caseSleepClean0 = (caseSleep1Clean0 + caseSleep2Clean0)/2;
caseSleepClean1 = (caseSleep1Clean1 + caseSleep2Clean1)/2;
caseSleepClean2 = (caseSleep1Clean2 + caseSleep2Clean2)/2;
controlWakeClean0 = (controlWake1Clean0 + controlWake2Clean0)/2;
controlWakeClean1 = (controlWake1Clean1 + controlWake2Clean1)/2;
controlWakeClean2 = (controlWake1Clean2 + controlWake2Clean2)/2;
controlSleepClean0 = (controlSleep1Clean0 + controlSleep2Clean0)/2;
controlSleepClean1 = (controlSleep1Clean1 + controlSleep2Clean1)/2;
controlSleepClean2 = (controlSleep1Clean2 + controlSleep2Clean2)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate box plots. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
groupIdx1 = [zeros(length(caseWake1Clean0),1); ones(length(controlWake1Clean0),1); ...
            ones(length(caseWake1Clean0),1)*2; ones(length(controlWake1Clean0),1)*3; ...
            ones(length(caseWake1Clean0),1)*4; ones(length(controlWake1Clean0),1)*5];

figure
subplot(1,2,1)
combinedData = [caseWakeClean0'; controlWakeClean0'; caseWakeClean1'; controlWakeClean1'; caseWakeClean2'; controlWakeClean2';];
boxplot(combinedData/19, groupIdx1, 'datalim', [0, 8.5], 'jitter', 0.1)
outlierTags = findobj(gca, 'tag', 'Outliers');
set(outlierTags, 'Marker', '.', 'markerEdgeColor', [1 0.5 0.5], 'markersize', 12)
ylabel('Canolty''s MI')
title('Awake')
set(gca, 'xtick', 1:6, 'xticklabel', {'Case0' 'Ctrl0' 'Case1' 'Ctrl1' 'Case2' 'Ctrl2'}, 'box', 'off')
ylim([0 9])

subplot(1,2,2)
combinedData = [caseSleepClean0'; controlSleepClean0'; caseSleepClean1'; controlSleepClean1'; caseSleepClean2'; controlSleepClean2';];
boxplot(combinedData/19, groupIdx1, 'datalim', [0, 8.5], 'jitter', 0.1)
outlierTags = findobj(gca, 'tag', 'Outliers');
set(outlierTags, 'Marker', '.', 'markerEdgeColor', [1 0.5 0.5], 'markersize', 12)
ylabel('Canolty''s MI')
title('Sleep')
set(gca, 'xtick', 1:6, 'xticklabel', {'Case0' 'Ctrl0' 'Case1' 'Ctrl1' 'Case2' 'Ctrl2'}, 'box', 'off')
ylim([0 9])

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)
set(gcf, 'position', [1 1 1986 1001])
print('/data/mobi/Hiroki/p4300_ageMatchedMis/boxplots', '-dsvg')

% Perform statistics.
[H,P,CI,STATS] = ttest2(caseWakeClean0', controlWakeClean0');

disp(sprintf('%.1f vs. %.1f', median(caseWakeClean0), median(controlWakeClean0)))
disp(sprintf('%.1f vs. %.1f', median(caseWakeClean1), median(controlWakeClean1)))
disp(sprintf('%.1f vs. %.1f', median(caseWakeClean2), median(controlWakeClean2)))

% 64.9 vs. 70.8
% 51.5 vs. 42.0
% 42.8 vs. 24.9

disp(sprintf('%.1f vs. %.1f', median(caseSleepClean0), median(controlSleepClean0)))
disp(sprintf('%.1f vs. %.1f', median(caseSleepClean1), median(controlSleepClean1)))
disp(sprintf('%.1f vs. %.1f', median(caseSleepClean2), median(controlSleepClean2)))

% 64.9 vs. 70.8
% 51.5 vs. 42.0
% 42.8 vs. 24.9

[P1,H1,STATS1] = ranksum(caseWakeClean0',  controlWakeClean0');
[P2,H2,STATS2] = ranksum(caseWakeClean1',  controlWakeClean1');
[P3,H3,STATS3] = ranksum(caseWakeClean2',  controlWakeClean2');
[P4,H4,STATS4] = ranksum(caseSleepClean0', controlSleepClean0');
[P5,H5,STATS5] = ranksum(caseSleepClean1', controlSleepClean1');
[P6,H6,STATS6] = ranksum(caseSleepClean2', controlSleepClean2');

% Compute Cohen's d. [0.2290    0.4574    1.1977    1.7224    1.5515   1.7364]
addpath('/data/projects/makoto/Tools')
output =   [cohens_d(caseWakeClean0, controlWakeClean0) ...
            cohens_d(caseWakeClean1, controlWakeClean1) ...
            cohens_d(caseWakeClean2, controlWakeClean2) ...
            cohens_d(caseSleepClean0, controlSleepClean0) ...
            cohens_d(caseSleepClean1, controlSleepClean1) ...
            cohens_d(caseSleepClean2, controlSleepClean2)];
        


%% Draw ROC curve.
groupIdx2 = [zeros(length(caseWake1Clean0),1); ones(length(controlWake1Clean0),1)];

[X1,Y1,T1,AUC1,OPTROCPT1] = perfcurve(groupIdx2, [caseWakeClean0'; controlWakeClean0'],   0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');
[X2,Y2,T2,AUC2,OPTROCPT2] = perfcurve(groupIdx2, [caseWakeClean1'; controlWakeClean1'],   0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');
[X3,Y3,T3,AUC3,OPTROCPT3] = perfcurve(groupIdx2, [caseWakeClean2'; controlWakeClean2'],   0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');
[X4,Y4,T4,AUC4,OPTROCPT4] = perfcurve(groupIdx2, [caseSleepClean0'; controlSleepClean0'], 0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');
[X5,Y5,T5,AUC5,OPTROCPT5] = perfcurve(groupIdx2, [caseSleepClean1'; controlSleepClean1'], 0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');
[X6,Y6,T6,AUC6,OPTROCPT6] = perfcurve(groupIdx2, [caseSleepClean2'; controlSleepClean2'], 0, 'xvals', 'all', 'nboot', 10000, 'BootType', 'normal');

disp(sprintf('%.3f vs. %.3f vs. %.3f', AUC1(1), AUC2(1), AUC3(1)))
disp(sprintf('%.3f vs. %.3f vs. %.3f', AUC4(1), AUC5(1), AUC6(1)))

% 0.431 vs. 0.665 vs. 0.845
% 1.000 vs. 0.987 vs. 0.984

figure
subplot(1,2,1)
plot(X1, Y1(:,1), 'linewidth', 2)
hold on
plot(X2, Y2(:,1), 'linewidth', 2)
plot(X3, Y3(:,1), 'linewidth', 2)
line([0 1], [0 1], 'color', [0 0 0])
ylim([0 1])
axis square
legend({'Cleaning Level 0', 'Cleaning Level 1', 'Cleaning Level 2'})
xlabel('1-Specificity')
ylabel('Sensitivity')
title('Awake')

subplot(1,2,2)
plot(X4, Y4(:,1), 'linewidth', 2)
hold on
plot(X5, Y5(:,1), 'linewidth', 2)
plot(X6, Y6(:,1), 'linewidth', 2)
line([0 1], [0 1], 'color', [0 0 0])
ylim([0 1])
axis square
legend({'Cleaning Level 0', 'Cleaning Level 1', 'Cleaning Level 2'})
xlabel('1-Specificity')
ylabel('Sensitivity')
title('Sleep')


set(findall(gcf, '-property', 'fontsize'), 'fontsize', 12)
set(gcf, 'position', [1 1 1986 1001])
print('/data/mobi/Hiroki/p4300_ageMatchedMis/rocCurves', '-dsvg')

% 
% figure
% plot(Y2)
% 
% figure
% plot(Y3)
% 
% load fisheriris
% x = meas(51:end,1:2);        % iris data, 2 classes and 2 features
% y = (1:100)'>50;             % versicolor=0, virginica=1
% b = glmfit(x,y,'binomial');  % logistic regression
% p = glmval(b,x,'logit');     % get fitted probabilities for scores
% 
% figure
% [X,Y] = perfcurve(species(51:end,:),p,'virginica');
% plot(X,Y)
% xlabel('False positive rate'); ylabel('True positive rate')
% title('ROC for classification by logistic regression')
% 
% % Obtain errors on TPR by vertical averaging
% [X,Y] = perfcurve(species(51:end,:),p,'virginica','nboot',10000,'xvals','all', 'BootType', 'normal');
% errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1)); % plot errors