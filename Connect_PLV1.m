function [CON_Matrix]=Connect_PLV(inte_data)
%信号间的相关
%input:
%global fs;1000hz
% global inte_data;它是去均值化后的15*10000*3的矩阵
%global Freq_limit;
%global Freq_num;
%output:
%global CON_Matrix;
% M=size(inte_data,1);%返回的是数组inte_data的行数为15
%h_w = waitbar(0,['正在计算网络连接强度，请稍候....']);进度条显示
%    
% for freq=1:1:Freq_num
% %     waitbar(freq/Freq_num);
%     low=Freq_limit{freq,1};
%     hig=Freq_limit{freq,2};
%     [filt_data] = eegfilt(inte_data,fs,low,hig,0,100);%数据的最小二乘发的滤波平滑处理
%     clc;
    A=data_psi(inte_data);%A是一个15*15的求了相位锁时值的矩阵
%     A(find(A==1))=0;%去1处理
    PLV = A;  %Phase synchronization index   
    %waitbar(freq/Freq_num);
% end
% CON_Matrix = [];
CON_Matrix = PLV;%15*15*3
%delete(h_w);

