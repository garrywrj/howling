load('HRIR_ANL_ha.mat')
hrtf=20*log10(abs(fft(hrir.impulseResponses(:,:,1))));
fs=48000;
N = fs/2/(length(hrtf)/2+1);
xi = 1:N:24000;
h=plot(xi,hrtf(1:129,:));
set(h,'LineWidth',2.2); 
xlabel('Ƶ�ʣ�Hz��','FontSize',12);
ylabel('����(dB)','FontSize',12);
title('HRTFƵ������');
legend('���','�Ҷ�');

% hrir=hrir.impulseResponses(:,:,1);
% h=plot(hrir);
% set(h,'LineWidth',2.2); 
% xlabel('������','FontSize',12);
% ylabel('����','FontSize',12);
% title('HRIRʱ������');
% legend('���','�Ҷ�');
% axis([0 256 -2 2]);