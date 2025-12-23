%% ==============================================
%% Network statistics + graph metrics (PLV-based)
%% ==============================================
clc; clear; close all;

%% ---------------- Parameters ----------------
name_band = ["\delta\","\theta\","\alpha\","\beta\","\gamma\"];
label_num = 1;          % frequency band index
bandName  = name_band(label_num);

nROI = 242;

%% ---------------- Paths ----------------
rootDir  = pwd; 
plvDir   = fullfile(rootDir,['a3_PLV', bandName]);
timeDir  = fullfile(rootDir,['a1_preprocessed', bandName]);

saveStat = fullfile(rootDir,['a5_plv_ttest', bandName]);
saveNet  = fullfile(rootDir,['a4_network_properties', bandName]);

mkdir(saveStat); mkdir(saveNet);

files = dir(fullfile(plvDir,'*.h5'));
nSub  = length(files);

%% ---------------- Load PLV data ----------------
R_plv   = zeros(nSub, nROI, nROI);
LOC_plv = zeros(nSub, nROI, nROI);

for s = 1:nSub
    fprintf('Loading subject %d / %d\n', s, nSub);

    % ---- load PLV ----
    plvFile = fullfile(plvDir, files(s).name);
    Data = hdf5read(plvFile, '/data');  % ROI ¡Á ROI ¡Á time

    % ---- load time marker ----
    name = files(s).name(1:end-3);
    load(fullfile(timeDir, name));      % variable: tt

    % ---- state definition ----
    Data_rest  = Data(:,:,116);                 % awake baseline
    Data_LOC   = Data(:,:,117+tt:end);          % post-LOC

    R_plv(s,:,:)   = mean(Data_rest,3);
    LOC_plv(s,:,:) = mean(Data_LOC,3);

    clear Data Data_rest Data_LOC
end

%% ---------------- Edge-wise statistics ----------------
fprintf('Running edge-wise statistics...\n');

nEdge = nROI * (nROI - 1) / 2;
alpha_corr = 0.001 / nEdge;

[hL,~,~,~] = ttest(R_plv, LOC_plv, alpha_corr, 'left');
[hR,~,~,~] = ttest(R_plv, LOC_plv, alpha_corr, 'right');

% Directional difference
B = squeeze(hR - hL);   % +1 LOC increase, -1 LOC decrease

% Keep upper triangle only
NetStat = zeros(nROI);
for i = 1:nROI
    for j = i+1:nROI
        NetStat(i,j) = B(i,j);
    end
end

save(fullfile(saveStat, ['PLV_edgeStats_' bandName '.mat']), 'NetStat','alpha_corr');

%% ---------------- Network visualization ----------------
NodenameFontSize = 9;
nodeMakerSize = 8;
roi_name={'FF','TT','PP','OO','Hipp','INS','CG','Tha','Cls'};
roi_color = [];
dong_NetPlot(NetStat, 2, 0, 0, roi_name, roi_color, NodenameFontSize, nodeMakerSize);

%% ---------------- Network properties ----------------
fprintf('Computing graph metrics...\n');
addpath(genpath('D:\MATLAB_toolbox\BCT_2015_01_25'));

Cluster_all  = zeros(nSub,1);
Charp_all    = zeros(nSub,1);
Cluster_LOC  = zeros(nSub,1);
Charp_LOC    = zeros(nSub,1);

for s = 1:nSub
    fprintf('Graph metrics: subject %d / %d\n', s, nSub);

    plvFile = fullfile(plvDir, files(s).name);
    Data = hdf5read(plvFile, '/data');  % ROI ¡Á ROI ¡Á time
    nWin = size(Data,3);

    cluster_sub  = zeros(nWin,1);
    charp_sub    = zeros(nWin,1);

    for t = 1:nWin
        W = Data(:,:,t);
        W(1:nROI+1:end) = 0;   % remove diagonal

        % ---- clustering coefficient ----
        Ci = clustering_coef_wu(W);
        cluster_sub(t) = mean(Ci);

        % ---- characteristic path length ----
        D = 1 - W;
        D(1:nROI+1:end) = 0;
        Dist = distance_wei(D);
        charp_sub(t) = sum(Dist(:)) / (nROI*(nROI-1));
    end

    % Average over awake vs LOC windows
    Cluster_all(s)  = mean(cluster_sub(1:116));
    Cluster_LOC(s)  = mean(cluster_sub(117+tt:end));
    Charp_all(s)    = mean(charp_sub(1:116));
    Charp_LOC(s)    = mean(charp_sub(117+tt:end));

    % Save individual metrics
    save(fullfile(saveNet, [files(s).name(1:end-3) '.mat']), ...
        'cluster_sub','charp_sub','Cluster_all','Cluster_LOC','Charp_all','Charp_LOC');

    clear Data cluster_sub charp_sub
end

%% ---------------- Network property statistics ----------------
fprintf('Running network property statistics...\n');

[~,p_clust,~,stats_clust] = ttest(Cluster_all, Cluster_LOC);
[~,p_charp,~,stats_charp] = ttest(Charp_all, Charp_LOC);

save(fullfile(saveNet, ['NetworkStats_' bandName '.mat']), ...
    'Cluster_all','Cluster_LOC','Charp_all','Charp_LOC', ...
    'p_clust','p_charp','stats_clust','stats_charp');

disp('===== All analysis finished =====');
