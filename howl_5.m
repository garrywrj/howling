clear all;
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
            for i=1:idx
                echoRow(i,:) = speakSignal(mod(i-1,buff)+1,1+2048*(mod(idx-i,buff)+1):2048*((mod(idx-i,buff)+1)+1));                
            end
            
        elseif idx>=buff
            speakSignal = [speakSignal(2:end,:);conv(farSpeech,roomImpulseResponse)'];
            
            for i=1:buff-1
                echoRow(i,:) = speakSignal(mod(i-1,buff)+1,1+2048*(mod(buff-1-i,buff)+1):2048*((mod(buff-1-i,buff)+1)+1));                
            end
        end
        
        farSpeechEcho = sum(echoRow,1);
       
end

%haimfigure;
plot(m);

% audiowrite('howl.wav',m,fs);
                
 

 


