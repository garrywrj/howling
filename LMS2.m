
[y1,fs]=audioread('C:\Users\garry\Desktop\howl.wav');
x = audioread('C:\Users\garry\Desktop\F_hecheng.wav');

T3=1000;                              %变步长自适应滤波器3滤波器阶数
w3=zeros(T3,1);                     %变步长自适应滤波器3抽头系数初始化
u3=0.0001;




for t=T3:length(y1)-T3                                      %变步长自适应滤波器3，变步长LMS过程
    fv3=y1(t:-1:t-T3+1);
    z3(t)=fv3'*w3;
    e3(t)=z3(t)-x(t);                                       %误差3
    u3=0.995*u3+0.0005*e3(t)*e3(t-1);                            %步长=误差3的函数
    w3=w3-2*u3*e3(t)*fv3;                                   %系数调整
end


plot(z3);
audiowrite('变步长lms2.wav',z3,fs);



