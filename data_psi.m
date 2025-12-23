function [PLV]=data_psi(data)
[M,N] = size(data);
if (M>N)
    data=data';   %一行为一导
    [M,N] = size(data);
end
PLV=zeros(size(data,1));  %15*15的0矩阵
for i=1:M-1     %m=15
     for j=(i+1):M
        PLV(i,j)=PhaseSI(data(i,:),data(j,:));
        PLV(j,i)=PLV(i,j);
     end 
end      %任意两行的相位锁时值得到一个15*15的矩阵，对角元为0
PLV=PLV+eye(size(data,1));%加上一个15*15的单位矩阵