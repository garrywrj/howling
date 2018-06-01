[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');


T=10000;                                           %自适应滤波器阶数
w1=zeros(T,1);                          %固定步长自适应滤波器1抽头初始化
u1=0.000005;                              %固定步长自适应滤波器1步长



for t=T:length(y1)-T                                        %固定步长自适应滤波器1，LMS过程
    fv1=y1(t:-1:t-T+1);                                     %滤波1
    z1(t)=fv1'*w1;                                          %输出1
    e1(t)=z1(t)-x(t);                                       %误差1
    w1=w1-2*u1*e1(t)*fv1;                                   %系数调整
end

plot(z1);
% audiowrite('lms.wav',z1,fs);
