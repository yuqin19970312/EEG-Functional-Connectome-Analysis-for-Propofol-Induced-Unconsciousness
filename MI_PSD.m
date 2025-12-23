function [PSD, SelectFilterFreq] = MI_PSD(x,FreqRange,fs)
%% 计算功率谱密度
% 输入：x：eeg数据，格式 time*channel
%           FreqRange：感兴趣的频率段范围，格式  for example，[5 30]
%           fs: 采样率
% 输出：PSD: 功率谱密度，格式 psd*channel
%           SelectFilterFreq：与psd相对应的频率点，格式 freq*channel
%
%        written by zhang rui in 2012.09.15
%

PSD_windowlength = fs;

Nfft=fs*2;
ChanNum=size(x,2);
%% 循环导联
for i=1:ChanNum
    [PSD(:,i) FilterFeq]=pwelch(x(:,i),PSD_windowlength,[],Nfft,fs);
end
[r,c,v]=find(PSD==inf);
if ~isempty(v)
    fprintf('row:%d col:%d is inf error\n',r,c);
end
FeqLow=find(FilterFeq>=FreqRange(1),1,'first');
FeqUp=find(FilterFeq<=FreqRange(2),1,'last');
FeqSeqRange=FeqLow:FeqUp;
PSD=PSD(FeqSeqRange,:);
SelectFilterFreq=FilterFeq(FeqSeqRange,:);

