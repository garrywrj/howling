[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');


T=10000;                                           %����Ӧ�˲�������
w1=zeros(T,1);                          %�̶���������Ӧ�˲���1��ͷ��ʼ��
u1=0.000005;                              %�̶���������Ӧ�˲���1����



for t=T:length(y1)-T                                        %�̶���������Ӧ�˲���1��LMS����
    fv1=y1(t:-1:t-T+1);                                     %�˲�1
    z1(t)=fv1'*w1;                                          %���1
    e1(t)=z1(t)-x(t);                                       %���1
    w1=w1-2*u1*e1(t)*fv1;                                   %ϵ������
end

plot(z1);
% audiowrite('lms.wav',z1,fs);
