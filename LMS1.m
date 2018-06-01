
[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');

T=10000;                                                    %自适应滤波器阶数
w2=zeros(T,1);                                              %变步长自适应滤波器2抽头初始化


for t=T:length(y1)-T                                        %变步长自适应滤波器2，变步长LMS过程
    fv2=y1(t:-1:t-T+1);
    z2(t)=fv2'*w2;
    e2(t)=z2(t)-x(t);                                       %误差2
    u2(t)=0.0001*(1-exp(-50*(e2(t)^2)));                    %步长=误差2的函数
    w2=w2-2*u2(t)*e2(t)*fv2;                                %系数调整
end


plot(z2);
audiowrite('变步长lms1.wav',z2,fs);