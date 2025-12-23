%% Source reconstruction and network construction (sLORETA + PLV)
clc; clear; close all;

%% ---------------- Parameters ----------------
fs      = 100;        % sampling rate after preprocessing
alpha   = 0.5;        % sLORETA regularization

bands = { ...
    'delta',[1 4]; ...
    'theta',[4 8]; ...
    'alpha',[8 13]; ...
    'beta',[13 30]; ...
    'gamma',[30 45]};

label_num = 3; % alpha band
bandName  = bands{label_num,1};
freqRange = bands{label_num,2};

%% ---------------- Paths ----------------
rootDir = pwd;  
dataDir = fullfile(rootDir,'a1_preprocessed',bandName);
saveSrc   = fullfile(rootDir,'a2_sLORETA_source', bandName);
saveROI   = fullfile(rootDir,'a2_ROI_timeseries', bandName);
savePLV   = fullfile(rootDir,'a3_PLV', bandName);
savePSD   = fullfile(rootDir,'a4_PSD', bandName);

mkdir(saveSrc); mkdir(saveROI); mkdir(savePLV); mkdir(savePSD);

%% ---------------- Load templates ----------------
load('lf_EGI128.mat'); % leadfield
ROI_Mark = textread('242-242-242-ROI-slor.txt');
ROI_Mark = ROI_Mark(1:242,:);

files = data_load(dataDir,'.mat');

%% ---------------- Loop over subjects ----------------
for sub = 1:length(files)

    fprintf('Subject %d / %d\n', sub, length(files));
    load(fullfile(dataDir, files{sub}));   % segments

    nSeg  = length(segments);
    nROI  = size(ROI_Mark,1);

    roi_ts = zeros(nROI, size(segments{1},2), nSeg);
    plv_mat  = zeros(nROI, nROI, nSeg);
    psd_mat  = [];

    for i = 1:nSeg

        %% Segment data
        eeg = segments{i};

        %% sLORETA source reconstruction
        src = wb_sloreta(eeg, lf, alpha);  % voxel ¡Á time

        %% ROI time series extraction
        roi_tc = zeros(nROI, size(src,2));
        for r = 1:nROI
            idx = ROI_Mark(r,:) > 0;
            roi_tc(r,:) = mean(src(idx,:), 1);
        end

        roi_ts(:,:,i) = roi_tc;

        %% PSD
        [psd_tmp, ~] = MI_PSD(roi_tc', freqRange, fs);
        psd_mat(:,:,i) = psd_tmp;

       %% PLV
        plv_mat(:,:,i) = Connect_PLV1(roi_tc);

       %% Save single-segment source
        h5file = fullfile(saveSrc, ...
            [files{sub}(1:end-4) '_seg' num2str(i) '.h5']);
        h5create(h5file, '/source', size(src));
        h5write(h5file, '/source', src);

        clear src roi_tc
    end

    %% Save subject-level data
    h5create(fullfile(saveROI,[files{sub}(1:end-4) '.h5']), '/roi', size(roi_ts));
    h5write (fullfile(saveROI,[files{sub}(1:end-4) '.h5']), '/roi', roi_ts);

    h5create(fullfile(savePLV,[files{sub}(1:end-4) '.h5']), '/plv', size(plv_mat));
    h5write (fullfile(savePLV,[files{sub}(1:end-4) '.h5']), '/plv', plv_mat);

    h5create(fullfile(savePSD,[files{sub}(1:end-4) '.h5']), '/psd', size(psd_mat));
    h5write (fullfile(savePSD,[files{sub}(1:end-4) '.h5']), '/psd', psd_mat);

    disp([bandName ' band finished for subject ' num2str(sub)]);
    clear segments roi_ts plv_mat psd_mat
end
