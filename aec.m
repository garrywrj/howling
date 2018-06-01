%% Acoustic Echo Cancellation (AEC)
% This example shows how to apply adaptive filters to acoustic echo
% cancellation (AEC).
% 
% Author(s): Scott C. Douglas
% Copyright 1999-2014 The MathWorks, Inc.

%% Introduction
% Acoustic echo cancellation is important for audio teleconferencing when
% simultaneous communication (or full-duplex transmission) of speech is
% necessary. In acoustic echo cancellation, a measured microphone signal
% |d(n)| contains two signals:
%
% * the near-end speech signal |v(n)|
% * the far-end echoed speech signal |dhat(n)|
%
% The goal is to remove the far-end echoed speech signal from the
% microphone signal so that only the near-end speech signal is transmitted.
% This example has some sound clips, so you might want to adjust your
% computer's volume now.

%% The Room Impulse Response
%
% First, we describe the acoustics of the loudspeaker-to-microphone signal
% path where the speakerphone is located. We can use a long finite impulse
% response filter to describe these characteristics. The following sequence
% of commands generates a random impulse response that is not unlike what a
% conference room would exhibit assuming a system sampling rate of |fs =
% 16000 Hz|.

fs = 16000;
M = fs/2 + 1;
frameSize = 8192;

[B,A] = cheby2(4,20,[0.1 0.7]);
IIR = dsp.IIRFilter('Numerator', [zeros(1,6) B], 'Denominator', A);

FVT = fvtool(IIR);  % Analyze the filter
FVT.Color = [1 1 1];

%%
H = step(IIR, ...
    (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
H = H/norm(H)*4;    % Room Impulse Response
firRoom = dsp.FIRFilter('Numerator', H');

fig = figure;
plot(0:1/fs:0.5, H);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Room Impulse Response');
fig.Color = [1 1 1];

%% The Near-End Speech Signal
%
% The teleconferencing system's user is typically located near the system's
% microphone. Here is what a male speech sounds like at the microphone.

v = audioread('C:\Users\garry\Desktop\test.wav');

AP              = dsp.AudioPlayer('SampleRate', fs);
nearSpeechSrc   = dsp.SignalSource('Signal',v,'SamplesPerFrame',frameSize);
nearSpeechScope = dsp.TimeScope('SampleRate', fs, ...
                    'TimeSpan', 4, ...
                    'YLimits', [-1.5 1.5], ...
                    'BufferLength', length(v), ...
                    'Title', 'Near-End Speech Signal', ...
                    'ShowGrid', true);

% Stream processing loop
while(~isDone(nearSpeechSrc))
    % Extract the speech samples from the input signal
    nearSpeech = step(nearSpeechSrc);
    % Send the speech samples to the output audio device
    step(AP, nearSpeech);
    % Plot the signal
    step(nearSpeechScope, nearSpeech);
end

%% The Far-End Speech Signal
%
% Now we describe the path of the far-end speech signal. A male voice
% travels out the loudspeaker, bounces around in the room, and then is
% picked up by the system's microphone. Let's listen to what his speech
% sounds like if it is picked up at the microphone without the near-end
% speech present.

x = 3*v;
farSpeechSrc    = dsp.SignalSource('Signal',x,'SamplesPerFrame',frameSize);
farSpeechSink   = dsp.SignalSink;
farSpeechScope  = dsp.TimeScope('SampleRate', fs, ...
                    'TimeSpan', 4, ...
                    'YLimits', [-5 5], ...
                    'BufferLength', length(x), ...
                    'Title', 'Far-End Speech Signal', ...
                    'ShowGrid', true);

% Stream processing loop
while(~isDone(farSpeechSrc))
    % Extract the speech samples from the input signal
    farSpeech = step(farSpeechSrc);
    % Add the room effect to the far-end speech signal
    farSpeechEcho = step(firRoom, farSpeech);
    % Send the speech samples to the output audio device
    step(AP, farSpeechEcho);
    % Plot the signal
    step(farSpeechScope, farSpeech);
    % Log the signal for further processing
    step(farSpeechSink, farSpeechEcho);
end

%% The Microphone Signal
%
% The signal at the microphone contains both the near-end speech and the
% far-end speech that has been echoed throughout the room. The goal of the
% acoustic echo canceler is to cancel out the far-end speech, such that
% only the near-end speech is transmitted back to the far-end listener.

reset(nearSpeechSrc);
farSpeechEchoSrc = dsp.SignalSource('Signal', farSpeechSink.Buffer, ...
                    'SamplesPerFrame', frameSize);
micSink         = dsp.SignalSink;
micScope        = dsp.TimeScope('SampleRate', fs,...
                    'TimeSpan', 4, ...
                    'YLimits', [-5 5], ...
                    'BufferLength', length(x), ...
                    'Title', 'Microphone Signal', ...
                    'ShowGrid', true);

% Stream processing loop
while(~isDone(farSpeechEchoSrc))
    % Microphone signal = echoed far-end + near-end + noise
    micSignal = step(farSpeechEchoSrc) + step(nearSpeechSrc) + ...
                    0.001*randn(frameSize,1);
    % Send the speech samples to the output audio device
    step(AP, micSignal);
    % Plot the signal
    step(micScope, micSignal);
    % Log the signal
    step(micSink, micSignal);
end

%% The Frequency-Domain Adaptive Filter (FDAF)
%
% The algorithm that we will use in this example is the *Frequency-Domain
% Adaptive Filter (FDAF)*. This algorithm is very useful when the impulse
% response of the system to be identified is long. The FDAF uses a fast
% convolution technique to compute the output signal and filter updates.
% This computation executes quickly in MATLAB(R). It also has improved
% convergence performance through frequency-bin step size normalization.
% We'll pick some initial parameters for the filter and see how well the
% far-end speech is cancelled in the error signal.

% Construct the Frequency-Domain Adaptive Filter
FDAF    = dsp.FrequencyDomainAdaptiveFilter('Length', 2048, ...
            'StepSize', 0.025, ...
            'InitialPower', 0.01, ...
            'AveragingFactor', 0.98, ...
            'Method', 'Unconstrained FDAF');

AECScope1   = dsp.TimeScope(4, fs, ...
                'LayoutDimensions', [4,1], ...
                'TimeSpan', 35, ...
                'BufferLength', length(x));

AECScope1.ActiveDisplay = 1;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1.5 1.5];
AECScope1.Title         = 'Near-End Speech Signal';

AECScope1.ActiveDisplay = 2;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1.5 1.5];
AECScope1.Title         = 'Microphone Signal';

AECScope1.ActiveDisplay = 3;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [-1.5 1.5];
AECScope1.Title         = 'Output of Acoustic Echo Canceller mu=0.025';

AECScope1.ActiveDisplay = 4;
AECScope1.ShowGrid      = true;
AECScope1.YLimits       = [0 50];
AECScope1.YLabel        = 'ERLE [dB]';
AECScope1.Title         = 'Echo Return Loss Enhancement mu=0.025';

% Near-end speech signal
release(nearSpeechSrc);
nearSpeechSrc.SamplesPerFrame = frameSize;

% Far-end speech signal
release(farSpeechSrc);
farSpeechSrc.SamplesPerFrame = frameSize;

% Far-end speech signal echoed by the room
release(farSpeechEchoSrc);
farSpeechEchoSrc.SamplesPerFrame = frameSize;

%% Echo Return Loss Enhancement (ERLE)
%
% Since we have access to both the near-end and far-end speech signals, we
% can compute the *echo return loss enhancement (ERLE)*, which is a
% smoothed measure of the amount (in dB) that the echo has been attenuated.
% From the plot, we see that we have achieved about a 35 dB ERLE at the end
% of the convergence period.

firERLE1 = dsp.FIRFilter('Numerator', ones(1,1024));
firERLE2 = clone(firERLE1);
setfilter(FVT,firERLE1);

micSrc = dsp.SignalSource('Signal', micSink.Buffer, ...
    'SamplesPerFrame', frameSize);

% Stream processing loop - adaptive filter step size = 0.025
while(~isDone(nearSpeechSrc))
    nearSpeech = step(nearSpeechSrc);
    farSpeech = step(farSpeechSrc);
    farSpeechEcho = step(farSpeechEchoSrc);
    micSignal = step(micSrc);
    % Apply FDAF
    [y,e] = step(FDAF, farSpeech, micSignal);
    % Send the speech samples to the output audio device
    step(AP, e);
    % Compute ERLE
    erle = step(firERLE1,(e-nearSpeech).^2)./ ...
        (step(firERLE2, farSpeechEcho.^2));
    erledB = -10*log10(erle);
    % Plot near-end, far-end, microphone, AEC output and ERLE
    step(AECScope1, nearSpeech, micSignal, e, erledB);
end

%% Effects of Different Step Size Values
%
% To get faster convergence, we can try using a larger step size value.
% However, this increase causes another effect, that is, the adaptive
% filter is "mis-adjusted" while the near-end speaker is talking. Listen
% to what happens when we choose a step size that is 60% larger than
% before.

% Change the step size value in FDAF
reset(FDAF);
FDAF.StepSize = 0.04;

AECScope2 = clone(AECScope1);
AECScope2.ActiveDisplay = 3;
AECScope2.Title = 'Output of Acoustic Echo Canceller mu=0.04';
AECScope2.ActiveDisplay = 4;
AECScope2.Title = 'Echo Return Loss Enhancement mu=0.04';

reset(nearSpeechSrc);
reset(farSpeechSrc);
reset(farSpeechEchoSrc);
reset(micSrc);
reset(firERLE1);
reset(firERLE2);

% Stream processing loop - adaptive filter step size = 0.04
while(~isDone(nearSpeechSrc))
    nearSpeech = step(nearSpeechSrc);
    farSpeech = step(farSpeechSrc);
    farSpeechEcho = step(farSpeechEchoSrc);
    micSignal = step(micSrc);
    % Apply FDAF
    [y,e] = step(FDAF, farSpeech, micSignal);
    % Send the speech samples to the output audio device
    step(AP, e);
    % Compute ERLE
    erle = step(firERLE1,(e-nearSpeech).^2)./ ...
        (step(firERLE2, farSpeechEcho.^2));
    erledB = -10*log10(erle);
    % Plot near-end, far-end, microphone, AEC output and ERLE
    step(AECScope2, nearSpeech, micSignal, e, erledB);
end

%% Echo Return Loss Enhancement Comparison
%
% With a larger step size, the ERLE performance is not as good due to the
% misadjustment introduced by the near-end speech. To deal with this
% performance difficulty, acoustic echo cancellers include a detection
% scheme to tell when near-end speech is present and lower the step size
% value over these periods. Without such detection schemes, the performance
% of the system with the larger step size is not as good as the former, as
% can be seen from the ERLE plots.

displayEndOfDemoMessage(mfilename)
