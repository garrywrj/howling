
[y1,fs]=audioread('C:\Users\garry\Desktop\howl.wav');
x = audioread('C:\Users\garry\Desktop\F_hecheng.wav');

T3=1000;                              %�䲽������Ӧ�˲���3�˲�������
w3=zeros(T3,1);                     %�䲽������Ӧ�˲���3��ͷϵ����ʼ��
u3=0.0001;




for t=T3:length(y1)-T3                                      %�䲽������Ӧ�˲���3���䲽��LMS����
    fv3=y1(t:-1:t-T3+1);
    z3(t)=fv3'*w3;
    e3(t)=z3(t)-x(t);                                       %���3
    u3=0.995*u3+0.0005*e3(t)*e3(t-1);                            %����=���3�ĺ���
    w3=w3-2*u3*e3(t)*fv3;                                   %ϵ������
end


plot(z3);
audiowrite('�䲽��lms2.wav',z3,fs);



