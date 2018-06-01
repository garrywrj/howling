clear;
%% The Room Impulse Response
% fs = 16000;
% M = fs/2 + 1;
rirLength = 5;
M = 2048*rirLength +1;
frameSize = 2048;

[B,A] = cheby2(4,20,0.7);%low pass filter
impulseResponseGenerator = dsp.IIRFilter('Numerator', [zeros(1,6) B], ...
    'Denominator', A);
roomImpulseResponse = impulseResponseGenerator( ...
        (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
roomImpulseResponse = roomImpulseResponse/norm(roomImpulseResponse);
% room = dsp.FIRFilter('Numerator', roomImpulseResponse');

%% The Near-End Speech Signal
v = audioread('F_hecheng.wav');

AudioInput = dsp.AudioFileReader(...
            'OutputDataType','double',...
            'Filename', 'F_hecheng.wav',...
            'PlayCount', 1);
fileInfo = info(AudioInput);
fs = fileInfo.SampleRate;
nearSpeechSrc   = dsp.SignalSource('Signal',v,'SamplesPerFrame',frameSize);

%% Gain
g = 4;

%% Initial the first-time echo by speaker
farSpeechEcho = 0;

%% Player
player          = audioDeviceWriter('SupportVariableSizeInput', true, ...
                                    'BufferSize', 512, 'SampleRate', fs);
                

%% write into a file
fileWrite = dsp.AudioFileWriter(...
	'test.wav',...
    'FileFormat','wav',...
	'SampleRate',fs);

%% lms setting
% sysorder = 5; % the number of filter orders

sysorder = 1000; % the number of filter orders
%% loop gain and loop phase 
% G = g;
% F = fft(roomImpulseResponse,frameSize);
% H = abs(G*F(1:frameSize/2+1));
% figure;
% plot(H);
% radP = phase(G*F(1:frameSize/2+1));
% figure;
% plot(1:length(radP),radP,'.');


%% Stream processing loop - adaptive filter step size = 0.025
speakerSignal = [];
idx = 0;
signalFar = [];
systemOutput = [];
rirRow = [];

while(~isDone(nearSpeechSrc))
        idx = idx + 1;
        
       % micSignal with effect of echo
        micSignal = nearSpeechSrc() + farSpeechEcho();
%         player(micSignal);

        % SpeakerSignal
        farSpeech = g * micSignal;        
%         player(farSpeech);

        % amplifier limit
        farSpeech((farSpeech)>1)=1;
        farSpeech((farSpeech)<-1)=-1;
%         step(fileWrite,farSpeech);

        % store the SpeakerSignal
        signalFar=[signalFar;farSpeech];
        
        % calculate roomImpluseResponse 
        nowRir = conv(farSpeech,roomImpulseResponse);
        if idx == 1
            rirRow = nowRir;
        else
            temp = bsxfun(@plus,rirRow(end-rirLength*frameSize+1:end,:),nowRir(1:rirLength*frameSize,:));
            rirRow = [rirRow(1:end-rirLength*frameSize);temp;nowRir(end-frameSize+1:end,:)];
        end
        
        %echo on next inputStream
        farSpeechEcho = rirRow(idx*frameSize+1:(idx+1)*frameSize);

        %lms processing        
        lmsInput = farSpeechEcho;               %input of lms,noise
        w=zeros(1,sysorder);                    %initiate the filter coefficient           
        err(1:sysorder) = micSignal(1:sysorder);  % The first N output assigned to its original value
        
        %calculate the output of lms
        for n = sysorder:frameSize            
            lmsOutput(n) = w * lmsInput(n-sysorder+1:1:n);  % lms filter output 
            err(n) = micSignal(n) - lmsOutput(n);                % the error of lms system, which also the output of system
%             if n<frameSize/2
%                  mu=0.001;
%              else
%                  mu=0.003;
%             end
            mu = 0.00002; 
            % update the coefficient
%             step = mu./(1+abs(err(n)).^2);
%             step = sqrt(step);
%             step = 2*mu;

%             u2(t)=0.0001*(1-exp(-50*(e2(t)^2))); 
            step = 0.0001*(1-exp(-50*abs(err(n)).^2));
            w = w + step*farSpeechEcho(n-sysorder+1:1:n)'*err(n); 

        end
            %the output of system
            oneOutput = err;        
            % store the system's output
            systemOutput = [systemOutput;oneOutput'];

        
       
end
figure;
plot(signalFar);

figure;
plot(systemOutput);
% player(systemOutput);
% step(fileWrite,systemOutput');
                
erledb=[];
Length = length(v);                
systemOutput=systemOutput(1:Length,:); 
m = v;
 for i=1:Length
    erledb(i) =10*log((m(i)^2)/(systemOutput(i)^2));
 end
erledb(isinf(erledb))=0;
erledb(isnan(erledb))=0;
ERLE = -mean(erledb);                
 
