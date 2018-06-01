
[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');


g=[0.000266;-0.002012;0.000897;-0.000504;0.000234;-0.000172;
   -0.000203;0.000293;-0.001131;-0.006891;-0.135781;-0.654303;
   -1.332799;-1.412939;-0.838869;-0.222403;0.065345;0.251710;
   0.296987;0.020201];                                       %�����ŵ�

h=fir2(200,[0,0.48,0.5,1],[1,1,0,0]);                        %���ŵ���Ӧ
h=h(:);
h=h.*exp(2*pi*i*(1:length(h))'/4);                           

h_dummy=zeros(size(h));
h_dummy((end+1)/2)=1;

k=0.3;                                                       %�������棬��Ҫ�������ı���Բ��Բ�������������

g=g(:);                                                      %����·��������
c=[0,0,0,0,1]';                                              %����ϵͳ�ڲ�����·��

xs1=zeros(size(c));
xs2=zeros(size(g));
xs3=zeros(size(h_dummy));
y1=zeros(size(x));
y2=zeros(size(x));

temp = 0;
f_shift=3;                                                    %��Ƶ3Hz



for j=1:length(x)
	xs1=[x(j)+temp;xs1(1:end-1)];
	y2(j)=k*(xs1'*c);

	xs3=[y2(j);xs3(1:end-1)];
	y2(j)=xs3'*h;                                            %�˲�������ź�Ƶ��������
	y2(j)=y2(j)*exp(2*pi*i*j/fs*f_shift);                    %�ƶ�Ƶ��
	y2(j)=real(y2(j));                                       %ȡʵ�������Ƶ�׸�����

	y2(j)=min(1,y2(j));                                      %�޷������y2����Ƶ����������Х���ź�
	y2(j)=max(-1,y2(j));
	xs2=[y2(j);xs2(1:end-1)];
	temp=xs2'*g;
end

plot(y2);
audiowrite('fshift.wav',y2,fs);






