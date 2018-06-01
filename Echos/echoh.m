tic;
[x,fs]=audioread('F_hecheng.wav');
x=x(1:fs*11);                                               %读取文件
g=[0.000266;-0.002012;0.000897;-0.000504;0.000234;-0.000172;
   -0.000203;0.000293;-0.001131;-0.006891;-0.135781;-0.654303;
   -1.332799;-1.412939;-0.838869;-0.222403;0.065345;0.251710;
   0.296987;0.020201];                                       %混响信道

h=fir2(200,[0,0.48,0.5,1],[1,1,0,0]);                        %功放的响应
h=h(:);
h=h.*exp(2*pi*i*(1:length(h))'/4);                           

h_dummy=zeros(size(h));
h_dummy((end+1)/2)=1;

k=0.3;                                                       %功放增益，重要参数，改变可以测试不出混响的最大功率

g=g(:);                                                      %反馈路径，混响
c=[0,0,0,0,1]';                                              %扩音系统内部传递路径

xs1=zeros(size(c));
xs2=zeros(size(g));
xs3=zeros(size(h_dummy));

y1=zeros(size(x));
y2=zeros(size(x));
temp=0;
 
for j=1:length(x)                                             %生成混响和啸叫，y1是生成信号
	xs1=[x(j)+temp;xs1(1:end-1)];                             %输入扩音器内部信号
	y1(j)=k*(xs1'*c);                                         %输出至扬声器信号 输入卷积内部路径，y1就是输出混响啸叫信号

	xs3=[y1(j);xs3(1:end-1)];
	y1(j)=xs3'*h_dummy;

	y1(j)=min(1,y1(j));                                       %幅度约束，出现啸叫就限定输出功率在-1和1之间
	y1(j)=max(-1,y1(j));
    
 	xs2=[y1(j);xs2(1:end-1)];                                 %扬声器输出信号
	temp=xs2'*g;                                              %扬声器输出信号卷积回音冲击形成混响
end

xs1=zeros(size(c));
xs2=zeros(size(g));
xs3=zeros(size(h));
temp=0;
f_shift=3;                                                    %移频3Hz

for j=1:length(x)
	xs1=[x(j)+temp;xs1(1:end-1)];
	y2(j)=k*(xs1'*c);

	xs3=[y2(j);xs3(1:end-1)];
	y2(j)=xs3'*h;                                            %滤波器获得信号频谱正半轴
	y2(j)=y2(j)*exp(2*pi*i*j/fs*f_shift);                    %移动频率
	y2(j)=real(y2(j));                                       %取实部，输出频谱负半轴

	y2(j)=min(1,y2(j));                                      %限幅，输出y2是移频处理后的抑制啸叫信号
	y2(j)=max(-1,y2(j));
	xs2=[y2(j);xs2(1:end-1)];
	temp=xs2'*g;
end


T=10000;                                                    %自适应滤波器阶数
u1=0.000005;                                                %固定步长自适应滤波器1步长
u3=0.00002;                                                 %变步长自适应滤波器2初始步长
w1=zeros(T,1);                                              %固定步长自适应滤波器1抽头初始化
w2=zeros(T,1);                                              %变步长自适应滤波器2抽头初始化


for t=T:length(y1)-T                                        %固定步长自适应滤波器1，LMS过程
    fv1=y1(t:-1:t-T+1);                                     %滤波1
    z1(t)=fv1'*w1;                                          %输出1
    e1(t)=z1(t)-x(t);                                       %误差1
    w1=w1-2*u1*e1(t)*fv1;                                   %系数调整
end


for t=T:length(y1)-T                                        %变步长自适应滤波器2，变步长LMS过程
    fv2=y1(t:-1:t-T+1);
    z2(t)=fv2'*w2;
    e2(t)=z2(t)-x(t);                                       %误差2
    u2(t)=0.0001*(1-exp(-50*(e2(t)^2)));     %步长=误差2的函数 beita=0.0001,alpha=50
    w2=w2-2*u2(t)*e2(t)*fv2;                                %系数调整
end
T3=1000;                                                    %变步长自适应滤波器3滤波器阶数
w3=zeros(T3,1);                                             %变步长自适应滤波器3抽头系数初始化
u3=0.0001;

for t=T3:length(y1)-T3                                      %变步长自适应滤波器3，变步长LMS过程
    fv3=y1(t:-1:t-T3+1);
    z3(t)=fv3'*w3;
    e3(t)=z3(t)-x(t);                                       %误差3
    u3=0.995*u3+0.0005*e3(t)*e3(t-1);                            %步长=误差3的函数
    w3=w3-2*u3*e3(t)*fv3;                                   %系数调整
end

figure(1),subplot(6,1,1),plot(x),title('输入信号'),ylabel('幅度');                             
          subplot(6,1,2),plot(y1),title(['啸叫混响信号，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,3),plot(y2),title(['移频法抑制啸叫后信号，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,4),plot(z1),title(['定步长LMS抑制啸叫后信号，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,5),plot(z2),title(['变步长LMS 1抑制啸叫后信号，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,6),plot(z3),title(['变步长LMS 2抑制啸叫后信号，增益k=',num2str(k)]),ylabel('幅度');
figure(2),subplot(6,1,1),plot(x),title('输入信号'),xlabel('采样点'),ylabel('幅度');
          subplot(6,1,2),plot(abs(y1-x)),title('输入信号'),ylabel('幅度');
          subplot(6,1,3),plot(abs(y2-x)),title(['移频法抑制啸叫误差，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,4),plot(abs(e1)),title(['定步长LMS抑制啸叫误差，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,5),plot(abs(e2)),title(['变步长LMS 1抑制啸叫误差，增益k=',num2str(k)]),ylabel('幅度');
          subplot(6,1,6),plot(abs(e3)),title(['变步长LMS 2抑制啸叫误差，增益k=',num2str(k)]),xlabel('采样点'),ylabel('幅度');
figure(3),stem(1:20,-g),title('混响冲击响应'),xlabel('采样点'),ylabel('幅度');
figure(4),freqz(h),title('功放滤波器幅频和相频响应');
figure(5),subplot(2,1,1),plot((-1024/2:1024/2-1)*fs/1024,abs(fftshift(fft(x,1024)))),title('输入信号频谱'),xlabel('幅度'),ylabel('频率');
          subplot(2,1,2),plot((-1024/2:1024/2-1)*fs/1024,abs(fftshift(fft(y1,1024)))),title('啸叫信号频谱'),xlabel('幅度'),ylabel('频率');

audiowrite('F_hechengecho.wav',y1,fs);
audiowrite('F_hechengecho_shift.wav',y2,fs);
audiowrite('F_hechengecho_lms.wav',z1,fs);
audiowrite('F_hechengecho_lms_ec.wav',z2,fs);
audiowrite('F_hechengecho_lms_rc.wav',z3,fs);
toc;

