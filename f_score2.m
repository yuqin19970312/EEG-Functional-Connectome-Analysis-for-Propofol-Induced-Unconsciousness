function F = f_score2(train_data,train_label)
 %train_data:训练集数据 样本×特征（电极）
 %train_label训练集标签 样本×1
%%----------计算特征的F值----------%%
Label = unique(train_label);
data1 =  train_data(train_label == Label(1,1),:); %第一类数据  9*44
%data2 =  train_data(train_label == Label(1,2),:); %第二类数据   根据数据不一样标签有可能不一样      21*44
data2 =  train_data(train_label == Label(2,1),:); %第二类数据    %根据数据不一样标签有可能不一样  21*44

%-----分别求出两组人特征的均值-----%
x1 = mean(data1); % 第一类特征的均值  44
x2 = mean(data2); % 第二类特征的均值  44
x_data = mean(train_data)  ; % 总体特征的均值 44

%-----F_score特征筛选-----% 
for i = 1:size(train_data,2) % 特征的个数 44
    for k = 1:size(train_data,1) % 样本的个数 30
            if k<=size(data1,1)  %9个第一类
                    AA(k) = (data1(k,i)-x1(i))^2;                                                                  %第1类人第k个人的第i个特征 - 第1类人第i个特征平均值    体现组内差异
                    sum1(i) = sum(AA);
              %  else if k<=(size(data1,1)+size(data2,1))  %30
                    else if k>size(data1,1) && k<=(size(data1,1)+size(data2,1))  %30
                    BB(k-size(data1,1)) = (data2((k-size(data1,1)),i)-x2(i))^2;                                    %第2类人第k个人的第i个特征 - 第2类人第i个特征平均值    体现组内差异
                    sum2(i) = sum(BB);
                    end
            end
    end
    F(i) = (((x1(i)-x_data(i))^2)+((x2(i)-x_data(i))^2))/((sum1(i)/(size(data1,1)-1))+(sum2(i)/(size(data2,1)-1))); %分子体现组间差异（枚举）
end