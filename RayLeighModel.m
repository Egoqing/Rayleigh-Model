%% 生成原始信息序列
clear
clc
simulation_point = 1000000;
R_b = 100000;%R_b=100kbps
F_s = 400000;%设定Frequency of Sampling = 400kbps
source_data = randsrc(1,simulation_point,[1,0]);%生成原始数据
%% 卷积编码
trellis = poly2trellis(9,[561 753]);
conv_data = convenc(source_data,trellis);
%% 交织
rows = 10000;cols = 10;%设定交织的深度与宽度
division = length(conv_data)/(rows*cols);%交织次数
interwine_data = zeros(1,2*simulation_point);
for i =1:division
    temp_data_1 = conv_data(1,(((i-1)*(rows*cols))+1):(i*(rows*cols)));
    temp_data_2 = matintrlv(temp_data_1,rows,cols);
    interwine_data(1:i*rows*cols) = horzcat(interwine_data(1:(i-1)*rows*cols),temp_data_2);
end
%% BPSK调制
send_data = interwine_data;
send_data(send_data ==0) = -1;%将01序列转换为1和-1发送
sample_data=zeros(1,4*simulation_point);
for i = 1:2*simulation_point
    sample_data((i-1)*2+1:i*2) = send_data(i);%要发送的序列
end
%% 单径Rayleigh信道
fd = 100; %最大多普勒频移
h=ones(1,2*simulation_point);
num=2*R_b/fd;%相干时间内编码器输出bit数
division2=2*simulation_point/(2*R_b/fd);
for i=1:division2
    temp_h = sqrt(1/2)*(randn(1,1)+j*randn(1,1));%复高斯分布
    h((i-1)*num+1:i*num)=abs(temp_h);%包络服从瑞利分布
end
for i = 1:2*simulation_point
    h_data((i-1)*2+1:i*2) = h(i);
end
EbN0 = 10;%取0~10dB 
N0 = 0.5*(F_s/(2*R_b))*10.^(-EbN0/10);%算出样本点噪声功率
noise = sqrt(N0).*randn(1,4*simulation_point);%算出样本点噪声功率
handle_data = sample_data.*h_data+noise; 
figure(1)
subplot(211);plot(sample_data(1:150));title('BPSK信号')
subplot(212);plot(handle_data(1:150));title('经过Rayleigh衰落后的信号')
%% 解调
receive_data=zeros(1,2*simulation_point);
for j = 1:(2*simulation_point)
    receive_data(j) = sum(handle_data(((j-1)*2+1):(j*2)));%匹配滤波？
end
%% 解交织
deinterwine_data = zeros(1,2*simulation_point);
for k = 1:division
    temp_data_1 = receive_data(1,(((k-1)*(rows*cols))+1):(k*(rows*cols)));
    temp_data_2 = matdeintrlv(temp_data_1,rows,cols);
    deinterwine_data(1:k*rows*cols) = horzcat(deinterwine_data(1:(k-1)*rows*cols),temp_data_2)';
end
%% 解码
Hdeinterwine_data=zeros(1,2*simulation_point);
Hdeinterwine_data(deinterwine_data>=0)=1;
Hdeinterwine_data(deinterwine_data<0)=0;
nsdec=3;%量化比特数
llim=min(min(deinterwine_data));
rlim=max(max(deinterwine_data));
partition=linspace(llim,rlim,2^nsdec+1);
partition(1)=[];partition(end)=[];
Sdeinterwine_data = quantiz(deinterwine_data,partition);
tbdepth=5;%A rate 1/2 code has a traceback depth of 5(ConstraintLength – 1).
Hdeconv_data=vitdec(Hdeinterwine_data,trellis,tbdepth,'trunc','hard');%硬判决
Sdeconv_data=vitdec(Sdeinterwine_data,trellis,tbdepth,'trunc','soft',3);%r软判决
figure(2)
subplot(311),stem(source_data(1:200));title('原始序列')
subplot(312),stem(Hdeconv_data(1:200));title('硬判决译码序列')
subplot(313),stem(Sdeconv_data(1:200));title('软判决译码序列')
%% 误码率计算
BERH = biterr(Hdeconv_data,source_data)/simulation_point;
BERS = biterr(Sdeconv_data,source_data)/simulation_point;
%% 错误点位置
indx=zeros(1,simulation_point);
indx(Hdeconv_data~=source_data)=1;
figure(3)
subplot(221);plot(indx(1:1000));title('前1/4的1000个点错误情况')
subplot(222);plot(indx(simulation_point/4+1:simulation_point/4+1000));title('前1/2的1000个点错误情况')
subplot(223);plot(indx(3*simulation_point/4+1:3*simulation_point/4+1000));title('前3/4的1000个点错误情况')
subplot(224);plot(indx(end-999:end));title('后1000个点错误情况')
SNR=10.^(EbN0/10);
BER_th=(1-sqrt(SNR/(1+SNR)))/2;