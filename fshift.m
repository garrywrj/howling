
[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');


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

temp = 0;
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

plot(y2);
audiowrite('fshift.wav',y2,fs);






