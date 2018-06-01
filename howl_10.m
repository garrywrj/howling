clear all;
%% The Room Impulse Response
fs = 16000;
% M = fs/2 + 1;
rirLength = 4;
M = 2048*rirLength +1;
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
farSpeechEchodum = 0;
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

%%
sysorder=5;%抽头数...

%% Stream processing loop - adaptive filter step size = 0.025
speakerSignal = [];
idx = 0;
echoRow = [];
signalFar = [];
signalFardum = [];
signalLsm = [];
rirRow = [];
rirRowdum = [];

% tst =0;
while(~isDone(nearSpeechSrc))
        idx = idx + 1;
        
       % micSignal with effect of echo
        micSignal = nearSpeechSrc() + farSpeechEcho();

%         
%   player(micSignal);
       
        % SpeakerSignal
        farSpeech = g * micSignal;
%         player(farSpeech);
        farSpeech((farSpeech)>1)=1;
        farSpeech((farSpeech)<-1)=-1;
        

        
        signalFar=[signalFar;farSpeech];

        % calculate roomImpluseResponse and plus the RIR to signal
        nowRir = conv(farSpeech,roomImpulseResponse);
        

        if idx == 1
            rirRow = nowRir;
        else
            temp = bsxfun(@plus,rirRow(end-rirLength*frameSize+1:end,:),nowRir(1:rirLength*frameSize,:));
            rirRow = [rirRow(1:end-rirLength*frameSize);temp;nowRir(end-frameSize+1:end,:)];
        end
        
        
  
        
    
        farSpeechEcho = rirRow(idx*frameSize+1:(idx+1)*frameSize);
       
        
%         
%         
%         ttt = nearSpeechSrc();
%         w=zeros(sysorder,1);%初始化
%         for n=sysorder:frameSize
%             y(1:sysorder)=micSignal(1:sysorder);  % 前五个值赋为原始值
%             y(n)=micSignal(n-sysorder+1:1:n)'*w;  %系统输出
%             e(n)=y(n)-ttt(n);  %系统误差
%             if n<200
%                  mu=0.0032;
%              else
%                  mu=0.0015;
%             end
%              step=mu./(1+abs(e(n)).^2);
%              w=w+step*micSignal(n:-1:n-sysorder+1)*e(n);%迭代方程
%         end
%         Y=-y;
%         player(Y');
%         signalLsm = [signalLsm;Y'];
%         
%        
end
% figure;
% plot(signalFar);
% figure;
% plot(signalLsm);



reset(nearSpeechSrc);
idx = 0;
f_shift = 4;
erle=[];
% for f_shift = 1:1:100
%     
while(~isDone(nearSpeechSrc))
        idx = idx + 1;
        
       % micSignal with effect of echo
        
        micSignaldum = nearSpeechSrc() + farSpeechEchodum();
  
%           player(micSignaldum);
        for j = 1 : frameSize
            micSignaldum(j) = micSignaldum(j)*exp(2*pi*1i*j/fs*f_shift);
            micSignaldum(j) = real(micSignaldum(j));
            micSignaldum(j) = min(1,micSignaldum(j));
            micSignaldum(j) = max(-1,micSignaldum(j)); 
        end
    


%         player(micSignal);
       
        % SpeakerSignal
        farSpeechdum = g * micSignaldum;


        farSpeechdum((farSpeechdum)>1)=1;
        farSpeechdum((farSpeechdum)<-1)=-1;
        
        signalFardum =[signalFardum;farSpeechdum];
        % calculate roomImpluseResponse and plus the RIR to signal
  
        nowRirdum = conv(farSpeechdum,roomImpulseResponse);
   
        if idx == 1
            rirRowdum = nowRirdum;
        else
            temp = bsxfun(@plus,rirRowdum(end-rirLength*frameSize+1:end,:),nowRirdum(1:rirLength*frameSize,:));
            rirRowdum = [rirRowdum(1:end-rirLength*frameSize);temp;nowRirdum(end-frameSize+1:end,:)];
        end
        
        
   
        
        farSpeechEchodum = rirRowdum(idx*frameSize+1:(idx+1)*frameSize);
                      
  
end
% answer = 0;
% for times = 1:1:180224
%     
%   answer =answer + (signalFardum(times)-signalFar(times))^2;
% end
% err(f_shift) = answer;
% reset(nearSpeechSrc);
% 
% end



% audiowrite('howl.wav',signalFar,fs);
% audiowrite('f_shift.wav',signalFardum,fs);

figure(1),subplot(3,1,1),plot(v);
          subplot(3,1,2),plot(signalFar);
          subplot(3,1,3),plot(signalFardum);
% figure;
% plot(signalLsm);


erledb=[];
Length = length(v);                
signalFardum=signalFardum(1:Length,:); 
m = v;
 for i=1:Length
    erledb(i) =10*log((m(i)^2)/(signalFardum(i)^2));
 end
erledb(isinf(erledb))=0;
erledb(isnan(erledb))=0;
ERLE = -mean(erledb);             








