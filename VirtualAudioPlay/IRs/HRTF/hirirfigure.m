load('HRIR_ANL_ha.mat')
hrtf=20*log10(abs(fft(hrir.impulseResponses(:,:,1))));
fs=48000;
N = fs/2/(length(hrtf)/2+1);
xi = 1:N:24000;
h=plot(xi,hrtf(1:129,:));
set(h,'LineWidth',2.2); 
xlabel('频率（Hz）','FontSize',12);
ylabel('幅度(dB)','FontSize',12);
title('HRTF频域曲线');
legend('左耳','右耳');

% hrir=hrir.impulseResponses(:,:,1);
% h=plot(hrir);
% set(h,'LineWidth',2.2); 
% xlabel('样点数','FontSize',12);
% ylabel('幅度','FontSize',12);
% title('HRIR时域曲线');
% legend('左耳','右耳');
% axis([0 256 -2 2]);