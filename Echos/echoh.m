tic;
[x,fs]=audioread('F_hecheng.wav');
x=x(1:fs*11);                                               %��ȡ�ļ�
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
temp=0;
 
for j=1:length(x)                                             %���ɻ����Х�У�y1�������ź�
	xs1=[x(j)+temp;xs1(1:end-1)];                             %�����������ڲ��ź�
	y1(j)=k*(xs1'*c);                                         %������������ź� �������ڲ�·����y1�����������Х���ź�

	xs3=[y1(j);xs3(1:end-1)];
	y1(j)=xs3'*h_dummy;

	y1(j)=min(1,y1(j));                                       %����Լ��������Х�о��޶����������-1��1֮��
	y1(j)=max(-1,y1(j));
    
 	xs2=[y1(j);xs2(1:end-1)];                                 %����������ź�
	temp=xs2'*g;                                              %����������źž����������γɻ���
end

xs1=zeros(size(c));
xs2=zeros(size(g));
xs3=zeros(size(h));
temp=0;
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


T=10000;                                                    %����Ӧ�˲�������
u1=0.000005;                                                %�̶���������Ӧ�˲���1����
u3=0.00002;                                                 %�䲽������Ӧ�˲���2��ʼ����
w1=zeros(T,1);                                              %�̶���������Ӧ�˲���1��ͷ��ʼ��
w2=zeros(T,1);                                              %�䲽������Ӧ�˲���2��ͷ��ʼ��


for t=T:length(y1)-T                                        %�̶���������Ӧ�˲���1��LMS����
    fv1=y1(t:-1:t-T+1);                                     %�˲�1
    z1(t)=fv1'*w1;                                          %���1
    e1(t)=z1(t)-x(t);                                       %���1
    w1=w1-2*u1*e1(t)*fv1;                                   %ϵ������
end


for t=T:length(y1)-T                                        %�䲽������Ӧ�˲���2���䲽��LMS����
    fv2=y1(t:-1:t-T+1);
    z2(t)=fv2'*w2;
    e2(t)=z2(t)-x(t);                                       %���2
    u2(t)=0.0001*(1-exp(-50*(e2(t)^2)));     %����=���2�ĺ��� beita=0.0001,alpha=50
    w2=w2-2*u2(t)*e2(t)*fv2;                                %ϵ������
end
T3=1000;                                                    %�䲽������Ӧ�˲���3�˲�������
w3=zeros(T3,1);                                             %�䲽������Ӧ�˲���3��ͷϵ����ʼ��
u3=0.0001;

for t=T3:length(y1)-T3                                      %�䲽������Ӧ�˲���3���䲽��LMS����
    fv3=y1(t:-1:t-T3+1);
    z3(t)=fv3'*w3;
    e3(t)=z3(t)-x(t);                                       %���3
    u3=0.995*u3+0.0005*e3(t)*e3(t-1);                            %����=���3�ĺ���
    w3=w3-2*u3*e3(t)*fv3;                                   %ϵ������
end

figure(1),subplot(6,1,1),plot(x),title('�����ź�'),ylabel('����');                             
          subplot(6,1,2),plot(y1),title(['Х�л����źţ�����k=',num2str(k)]),ylabel('����');
          subplot(6,1,3),plot(y2),title(['��Ƶ������Х�к��źţ�����k=',num2str(k)]),ylabel('����');
          subplot(6,1,4),plot(z1),title(['������LMS����Х�к��źţ�����k=',num2str(k)]),ylabel('����');
          subplot(6,1,5),plot(z2),title(['�䲽��LMS 1����Х�к��źţ�����k=',num2str(k)]),ylabel('����');
          subplot(6,1,6),plot(z3),title(['�䲽��LMS 2����Х�к��źţ�����k=',num2str(k)]),ylabel('����');
figure(2),subplot(6,1,1),plot(x),title('�����ź�'),xlabel('������'),ylabel('����');
          subplot(6,1,2),plot(abs(y1-x)),title('�����ź�'),ylabel('����');
          subplot(6,1,3),plot(abs(y2-x)),title(['��Ƶ������Х��������k=',num2str(k)]),ylabel('����');
          subplot(6,1,4),plot(abs(e1)),title(['������LMS����Х��������k=',num2str(k)]),ylabel('����');
          subplot(6,1,5),plot(abs(e2)),title(['�䲽��LMS 1����Х��������k=',num2str(k)]),ylabel('����');
          subplot(6,1,6),plot(abs(e3)),title(['�䲽��LMS 2����Х��������k=',num2str(k)]),xlabel('������'),ylabel('����');
figure(3),stem(1:20,-g),title('��������Ӧ'),xlabel('������'),ylabel('����');
figure(4),freqz(h),title('�����˲�����Ƶ����Ƶ��Ӧ');
figure(5),subplot(2,1,1),plot((-1024/2:1024/2-1)*fs/1024,abs(fftshift(fft(x,1024)))),title('�����ź�Ƶ��'),xlabel('����'),ylabel('Ƶ��');
          subplot(2,1,2),plot((-1024/2:1024/2-1)*fs/1024,abs(fftshift(fft(y1,1024)))),title('Х���ź�Ƶ��'),xlabel('����'),ylabel('Ƶ��');

audiowrite('F_hechengecho.wav',y1,fs);
audiowrite('F_hechengecho_shift.wav',y2,fs);
audiowrite('F_hechengecho_lms.wav',z1,fs);
audiowrite('F_hechengecho_lms_ec.wav',z2,fs);
audiowrite('F_hechengecho_lms_rc.wav',z3,fs);
toc;

