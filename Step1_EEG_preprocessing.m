%% EEG preprocessing for propofol-induced unconsciousness

clear; clc;clear all;

%% ---------------- parameters ----------------
fs_raw = 500;        % original sampling rate
fs = 100;            % downsampled rate
notchBand = [49 51]; % 50 Hz notch
fullBand = [1 45];   % broadband
bands = { ...
    'delta',[1 4]; ...
    'theta',[4 8]; ...
    'alpha',[8 13]; ...
    'beta',[13 30]; ...
    'gamma',[30 45]};

winSec = 5;          % window length (s)
overlap = 0.8;       % 80% overlap

%% ---------------- Paths ----------------
rawDir  = 'DATA_ROOT/a0_raw_mff/'; 
saveDir = 'DATA_ROOT/a1_preprocessed/';

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

subFiles = dir(fullfile(rawDir,'*.mff'));

%% ---------------- Loop over subjects ----------------
for s = 1:length(subFiles)

    fprintf('Processing subject %d / %d\n',s,length(subFiles));

    % Load raw EEG
    EEG = mff_import(fullfile(rawDir,subFiles(s).name));

    % 50 Hz notch filter
    EEG = pop_eegfiltnew(EEG, notchBand(1), notchBand(2), [], 1);

    % 1¨C45 Hz bandpass filter
    EEG = pop_eegfiltnew(EEG, fullBand(1), fullBand(2));

    % Load manually identified bad channels
    load(fullfile('BAD_CHANNELS', [subFiles(s).name '_badchans.mat']));

    % Exclude subject if >10 bad channels
    if length(badChans) > 10
        fprintf('Subject excluded (>10 bad channels)\n');
        continue
    end
    
    % REST re-referencing
    addpath(genpath('REST_v1.2_20200818'));
    EEG = rest_refer(EEG);

    % Interpolate bad channels
    EEG = pop_interp(EEG, badChans, 'spherical');

    % ICA (Infomax)
    EEG = pop_runica(EEG, 'extended', 1, 'stop', 1e-7);

    % IC classification and rejection
    EEG = pop_iclabel(EEG,'default');
    EEG = processMARA(EEG);

    % Reject artifact components (combined ICLabel + MARA + visual check)
    EEG = pop_icflag(EEG, ...
        [NaN NaN; ...  % Brain
         0.9 1;  ...  % Muscle
         NaN NaN; ...  % Eye
         NaN NaN; ...  % Heart
         NaN NaN; ...  % Line noise
         NaN NaN; ...  % Channel noise
         NaN NaN]);    % Other

    % Down-sample to 100 Hz
    EEG = pop_resample(EEG, fs);

    % Extract anesthesia-centered 4-min segment
    % Event extraction
    point = [EEG.event.latency];
    mark  = {EEG.event.type};

    s_inj = []; n_loss = [];
    for k = 1:length(mark)
        if strcmp(mark{k},'SEIZ')   % anesthetic injection
            s_inj = point(k);
        elseif strcmp(mark{k},'SPIK') % loss of consciousness
            n_loss = point(k);
        end
    end

    tt = ceil((n_loss - s_inj)/fs);
    EEG.data  = EEG.data(:, s_inj-120*fs : s_inj+120*fs);
    EEG.pnts  = size(EEG.data,2);
    EEG.times = EEG.times(1:EEG.pnts);

    %% Band-specific filtering and sliding-window segmentation
    for b = 1:size(bands,1)

        bandName = bands{b,1};
        bandRange = bands{b,2};

        EEGb = pop_eegfiltnew(EEG, bandRange(1), bandRange(2), 200);

        win = winSec * fs;
        step = win * (1-overlap);
        numWin = floor((size(EEGb.data,2)-win)/step);

        segments = cell(numWin,1);
        for w = 1:numWin
            idx = round((w-1)*step + (1:win));
            segments{w} = EEGb.data(:,idx);
        end

        save(fullfile(saveDir, ...
            [subFiles(s).name '_' bandName '.mat']), ...
            'segments','tt','-v7.3');

        clear EEGb segments
    end

    clear EEG
end
