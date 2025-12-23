%% ===============================================
%% SVM classification based on PLV features
%% ===============================================
clear; clc; close all;

%% ---------------- Parameters ----------------
name_band   = ["\delta\","\theta\","\alpha\","\beta\","\gamma\"];

%% ---------------- Initialize figure ----------------
figure; hold on;

%% ---------------- Loop over frequency bands ----------------
for b = 1:length(name_band)
    
    fprintf('Processing %s band...\n', namee(b));
    
    % ---- Load PLV-based FC features ----
    RR_data  = [];  LOC_data  = [];
    load(fullfile([name_band(b),'Rest.mat'])); % variable CC
    load(fullfile( [name_band(b),'LOC.mat']));  % variable DD
    RR_data  = [RR_data; CC];
    LOC_data = [LOC_data; DD];
    
    % ---- Labels ----
    nRR  = size(RR_data,2);
    nLOC = size(LOC_data,2);
    label = [ones(nRR,1); -ones(nLOC,1)];  % 1 = Rest, -1 = LOC
    
    % ---- Combine features ----
    train_feature = [RR_data LOC_data]';  % features ¡Á samples
    train_feature = train_feature(:,idx); % reorder columns
    
    %% ---------------- LOOCV ----------------
    nSamples = size(train_feature,1);
    pre_label = zeros(nSamples,1);
    score     = zeros(nSamples,2);
    
    for i = 1:nSamples
        % Training data
        Train_Feature = train_feature; Train_Feature(i,:) = [];
        Label_Train   = label;        Label_Train(i) = [];
        Test_Feature  = train_feature(i,:);
        Label_Test    = label(i);
        
        % ---- F-score feature selection (top 25%) ----
        F = f_score2(Train_Feature, Label_Train);
        [~, IX_sort] = sort(F,'descend');
        K = floor(length(F)*0.25);
        Train_Feature_sel = Train_Feature(:, IX_sort(1:K));
        Test_Feature_sel  = Test_Feature(:, IX_sort(1:K));
        
        % ---- SVM parameter grid search on training set ----
        C_list = [0.1 1 10];
        gamma_list = [0.1 0.5 1];
        bestAcc = 0; bestC = 1; bestGamma = 0.5;
        for C = C_list
            for gamma = gamma_list
                svmModel = fitcsvm(Train_Feature_sel, Label_Train, ...
                                   'KernelFunction','rbf', ...
                                   'BoxConstraint',C, ...
                                   'KernelScale',gamma, ...
                                   'Standardize',true, ...
                                   'CrossVal','on','KFold',5);  % 5-fold CV in training
                acc = 1 - kfoldLoss(svmModel);
                if acc > bestAcc
                    bestAcc = acc;
                    bestC = C; bestGamma = gamma;
                end
            end
        end
        
        % ---- Train final SVM with best parameters ----
        svmModel = fitcsvm(Train_Feature_sel, Label_Train, ...
                           'KernelFunction','rbf', ...
                           'BoxConstraint',bestC, ...
                           'KernelScale',bestGamma, ...
                           'Standardize',true);
        
        % ---- Predict test sample ----
        [pre_label(i), score(i,:)] = predict(svmModel, Test_Feature_sel);
    end
    
    %% ---------------- Performance metrics ----------------
    Acc = 100*mean(pre_label==label);
    Sen = 100*sum(pre_label(nRR+1:end)==-1)/nLOC; % Sensitivity: LOC
    Spe = 100*sum(pre_label(1:nRR)==1)/nRR;       % Specificity: Rest
    
    fprintf('%s band: Accuracy=%.2f%%, Sensitivity=%.2f%%, Specificity=%.2f%%\n', ...
        namee(b), Acc, Sen, Spe);
    
    %% ---------------- ROC ----------------
    [X_roc, Y_roc, ~, AUC(b)] = perfcurve(label, score(:,2), 1);
    plot(X_roc, Y_roc, 'Color', color(b,:), 'LineWidth', 1);
end

%% ---------------- Plot ROC ----------------
plot(0:0.01:1, 0:0.01:1, 'k--', 'LineWidth', 1);
legend(...
    ['Delta, AUC=',num2str(AUC(1))], ...
    ['Theta, AUC=',num2str(AUC(2))], ...
    ['Alpha, AUC=',num2str(AUC(3))], ...
    ['Beta, AUC=',num2str(AUC(4))], ...
    ['Gamma, AUC=',num2str(AUC(5))], ...
    'Location','southeast');

xlabel('1-Specificity'); ylabel('Sensitivity');
title('SVM classification ROC');
xticks(0:0.1:1); yticks(0:0.1:1); grid on; axis([0 1 0 1]);
