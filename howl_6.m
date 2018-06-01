clear all;
tic;
%% The Room Impulse Response
fs = 16000;
% M = fs/2 + 1;
M = 2048*4 +1;
frameSize = 2048;

[B,A] = cheby2(4,20,[0.1 0.7]);
impulseResponseGenerator = dsp.IIRFilter('Numerator', [zeros(1,6) B], ...
    'Denominator', A);

% FVT = fvtool(impulseResponseGenerator);  % Analyze the filter
% FVT.Color = [1 1 1];

%%
roomImpulseResponse = impulseResponseGenerator( ...
        (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
roomImpulseResponse = roomImpulseResponse/norm(roomImpulseResponse)*4;
room = dsp.FIRFilter('Numerator', roomImpulseResponse');

% fig = figure;
% plot(0:1/fs:0.5, roomImpulseResponse);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Room Impulse Response');
% fig.Color = [1 1 1];

%% The Near-End Speech Signal
v = audioread('F_hecheng.wav');

AudioInput = dsp.AudioFileReader(...
            'OutputDataType','double',...
            'Filename', 'F_hecheng.wav',...
            'PlayCount', 1);
fileInfo = info(AudioInput);
fs = fileInfo.SampleRate;
nearSpeechSrc   = dsp.SignalSource('Signal',v,'SamplesPerFrame',frameSize);
% nearSpeechScope = dsp.TimeScope('SampleRate', fs, ...
%                     'TimeSpan', 15, 'TimeSpanOverrunAction', 'Scroll', ...
%                     'YLimits', [-1.5 1.5], ...
%                     'BufferLength', length(v), ...
%                     'Title', 'Near-End Speech Signal', ...
%                     'ShowGrid', true);

%% Gain
g = 3;

%% The Far-End Speech Signal
farSpeechSink   = dsp.SignalSink;
% farSpeechScope  = dsp.TimeScope('SampleRate', fs, ...
%                     'TimeSpan', 35, 'TimeSpanOverrunAction', 'Scroll', ...
%                     'YLimits', [-1*g 1*g], ...
%                     'BufferLength', length(v), ...
%                     'Title', 'Far-End Speech Signal', ...
%                     'ShowGrid', true);

%% Initial the first-time echo by speaker
farSpeechEcho = 0;

%% Player
player          = audioDeviceWriter('SupportVariableSizeInput', true, ...
                                    'BufferSize', 512, 'SampleRate', fs);

                
%% The Frequency-Domain Adaptive Filter (FDAF)               
% Construct the Frequency-Domain Adaptive Filter
echoCanceller    = dsp.FrequencyDomainAdaptiveFilter('Length', 2048, ...
                    'StepSize', 0.025, ...
                    'InitialPower', 0.01, ...
                    'AveragingFactor', 0.98, ...
                    'Method', 'Unconstrained FDAF');  
                
                

               
%% time delay convert to samples                  
%	
timeDelayOfPowerAmplifier = 0.000;  %1ms 
IntervalOfSamples = 1/fs;
samplesOfPowerAmplifier = ceil(timeDelayOfPowerAmplifier/IntervalOfSamples);  %1ms 
sPA = samplesOfPowerAmplifier;


%% Stream processing loop - adaptive filter step size = 0.025
speakerSignal = [];
idx = 0;
buff = 5;
frameNum = ceil(length(v)/frameSize);
echoRow=[];
m=[];
while(~isDone(nearSpeechSrc))
        idx = idx + 1;
        
       % micSignal with effect of echo
        micSignal = nearSpeechSrc() + farSpeechEcho()';
%         player(micSignal);


        % Power limiter
%         if ((micSignal'*micSignal)> 10)
%             g = g/sqrt(micSignal'*micSignal);
%         end
        
        % SpeakerSignal
        farSpeech = g * micSignal;
        player(farSpeech);
        m=[m;farSpeech];
        
        % store the output signal
        if idx<buff
            speakSignal(mod(idx-1,buff)+1,:) = conv(farSpeech,roomImpulseResponse);            
            for j=1:idx
                echoRow(j,:) = speakSignal(mod(j-1,buff)+1,1+2048*(mod(idx-j,buff)+1):2048*((mod(idx-j,buff)+1)+1));                
            end
            
        elseif idx>=buff
            speakSignal = [speakSignal(2:end,:);conv(farSpeech,roomImpulseResponse)'];
            
            for j=1:buff-1
                echoRow(j,:) = speakSignal(mod(j-1,buff)+1,1+2048*(mod(buff-1-j,buff)+1):2048*((mod(buff-1-j,buff)+1)+1));                
            end
        end
        
        farSpeechEcho = sum(echoRow,1);
       
end
 g=[0.000266;-0.002012;0.000897;-0.000504;0.000234;-0.000172;
   -0.000203;0.000293;-0.001131;-0.006891;-0.135781;-0.654303;
   -1.332799;-1.412939;-0.838869;-0.222403;0.065345;0.251710;
   0.296987;0.020201];                                     
h=fir2(200,[0,0.48,0.5,1],[1,1,0,0]);                      
h=h(:);
h=h.*exp(2*pi*i*(1:length(h))'/4);                           
h_dummy=zeros(size(h));
h_dummy((end+1)/2)=1;
k=0.3;                                                      
g=g(:);                                                      
c=[0,0,0,0,1]';                                             
xs1=zeros(size(c));
xs2=zeros(size(g));
xs3=zeros(size(h_dummy));
y2=zeros(size(v));
%haimfigure; 


                
% 




temp = 0;
f_shift=3;                                                    %移频3Hz



for j=1:length(v)
	xs1=[v(j)+temp;xs1(1:end-1)];
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


[y1,fs]=audioread('howl.wav');


T=10000;                                           %自适应滤波器阶数
w1=zeros(T,1);                          %固定步长自适应滤波器1抽头初始化
u1=0.000005;                              %固定步长自适应滤波器1步长



for t=T:length(y1)-T                                        %固定步长自适应滤波器1，LMS过程
    fv1=y1(t:-1:t-T+1);                                     %滤波1
    z1(t)=fv1'*w1;                                          %输出1
    e1(t)=z1(t)-v(t);                                       %误差1
    w1=w1-2*u1*e1(t)*fv1;                                   %系数调整
end




for j=1:length(m)                                           
 	m(j)=min(1,m(j));                                       %幅度约束，出现啸叫就限定输出功率在-1和1之间
 	m(j)=max(-1,m(j));
end

figure(1),subplot(4,1,1),plot(v),title('输入信号'),ylabel('幅度');                             
          subplot(4,1,2),plot(m),title('啸叫信号'),ylabel('幅度');
          subplot(4,1,3),plot(y2),title('移频法抑制啸叫后信号'),ylabel('幅度');
          subplot(4,1,4),plot(z1),title('LMS抑制啸叫后信号'),ylabel('幅度');
audiowrite('howl.wav',m,fs);
audiowrite('fshift.wav',y2,fs);
audiowrite('lms.wav',z1,fs);

toc;

