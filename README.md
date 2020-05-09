电子科技大学-信通学院-移动通信课程作业
# 作业要求
单径瑞利信道；最大多普勒频移100Hz。<br>
信源比特速率R_b:100Kbps。<br>
卷积编码：码率1/2；八进制生成多项式(561，753)。<br>
交织：行列交织，深度100<br>
仿真点数1,000,000,000。<br>
译码方式硬判决译码与软判决译码(3bits)。<br>
Q1:AWGN信道下硬判决与软判决误码率曲线仿真。<br>Q2：单径瑞利信道下(最大多普勒频移100Hz)有无交织，不同译码条件下的误码率曲线仿真。<br>Q3：不同深度下(自拟)的误码率曲线仿真<br>
# 实现
## 前言
在无编码无交织条件下，AWGN信道与瑞利信道理论误码率为：<br>
![理论误码率曲线](https://github.com/Egoqing/Rayleigh-Model/blob/master/result/theoretical.jpg)<br>
## Q1
卷积码可以用MATLAB函数poly2trellis与convenc实现。硬判决与软判决维特比译码可以用poly2trellis函数实现。<br>
交织用matlab自带的matintrlv函数实现。<br>
利用Eb/NO计算噪声功率N0公式(实信号):N0 = 0.5*(F_s/R_b)*10.^(-Eb/N0/10)。<br>
噪声序列生成：noise= sqrt(N0).*randn(1,4*simulation_point)。<br>
![AWGN信道下仿真误码率曲线](https://github.com/Egoqing/Rayleigh-Model/blob/master/result/q1.jpg)<br>
