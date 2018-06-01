
[y1,fs]=audioread('howl.wav');
x = audioread('F_hecheng.wav');

T=10000;                                                    %����Ӧ�˲�������
w2=zeros(T,1);                                              %�䲽������Ӧ�˲���2��ͷ��ʼ��


for t=T:length(y1)-T                                        %�䲽������Ӧ�˲���2���䲽��LMS����
    fv2=y1(t:-1:t-T+1);
    z2(t)=fv2'*w2;
    e2(t)=z2(t)-x(t);                                       %���2
    u2(t)=0.0001*(1-exp(-50*(e2(t)^2)));                    %����=���2�ĺ���
    w2=w2-2*u2(t)*e2(t)*fv2;                                %ϵ������
end


plot(z2);
audiowrite('�䲽��lms1.wav',z2,fs);