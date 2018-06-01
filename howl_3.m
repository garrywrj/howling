clear all;
%% The Room Impulse Response
fs = 16000;
M = fs/2 + 1;
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
g = 1.2;

%% The Far-End Speech Signal
farSpeechSink   = dsp.SignalSink;
% farSpeechScope  = dsp.TimeScope('SampleRate', fs, ...
%                     'TimeSpan', 35, 'TimeSpanOverrunAction', 'Scroll', ...
%                     'YLimits', [-1*y 1*y], ...
%                     'BufferLength', length(v), ...
%                     'Title', 'Far-End Speech Signal', ...
%                     'ShowGrid', true);

%% Initial the first-time echo by speaker
farSpeechEcho = 0;

%% The micSignal

%% Player
player          = audioDeviceWriter('SupportVariableSizeInput', true, ...
                                    'BufferSize', 512, 'SampleRate', fs);
%% Buffer
farSpeechEchoBuff = [];
farSpeechBuff = [];

%% Stream processing loop for the first time 
% check we get howl or not
% no need to process the howl 

while(~isDone(nearSpeechSrc))
        % micSignal with effect of echo
        micSignal = nearSpeechSrc() + farSpeechEcho;
        % micSignal with no-effect of echo
%         micSignal = nearSpeechSrc();
        % SpeakerSignal
        farSpeech = g * micSignal;
        % Send the speech samples to the output audio device
        player(farSpeech);
        % Plot the  signal
%         farSpeechScope(farSpeech);
        % Add the room effect to the far-end speech signal
        farSpeechEcho = room(farSpeech);

        % Log the signal for further processing
%         farSpeechSink(farSpeechEcho);
end
                
%% The Frequency-Domain Adaptive Filter (FDAF)               
% Construct the Frequency-Domain Adaptive Filter
% echoCanceller    = dsp.FrequencyDomainAdaptiveFilter('Length', 2048, ...
%                     'StepSize', 0.025, ...
%                     'InitialPower', 0.01, ...
%                     'AveragingFactor', 0.98, ...
%                     'Method', 'Unconstrained FDAF');
                   
%% Stream processing loop - adaptive filter step size = 0.025
% while(~isDone(nearSpeechSrc))
%         % micSignal with effect of echo
%         micSignal = nearSpeechSrc() + farSpeechEcho;
%         % SpeakerSignal
%         farSpeech = g * micSignal;
%         % Add the room effect to the far-end speech signal
%         farSpeechEcho = room(farSpeech);        
%         % Apply FDAF
%         [y,e] = echoCanceller(farSpeech, micSignal);
%         % Send the speech samples to the output audio device
%         player(e);
%     
%         % Compute ERLE
% %     erle = diffAverager((e-nearSpeech).^2)./ farEchoAverager(farSpeechEcho.^2);
% %     erledB = -10*log10(erle);
%     % Plot near-end, far-end, microphone, AEC output and ERLE
% %     AECScope1(nearSpeech, micSignal, e, erledB);
% end                
%                 
                
 

 


