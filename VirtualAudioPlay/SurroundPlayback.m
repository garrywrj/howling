clc
clear all

% addpath('..\HOAToolbox');

% head tracker on/off
HT = 1;
R = 100; % Hz
Tim = 1./R; % in seconds

% individual HRTF
indHrtf = true;

if indHrtf == 0
    subjIni = 'ha';
else
    % Get listener's name
    listenerId = inputdlg('Please enter your name.', 'Comment',[1,40]) ;
    subjIni = char(listenerId) ;
end

% HRTF/DTF mode
dtf = false ;

% Audio parameters: buffer size, Head Tracker on/off ...
buffLen = 512 ;
buffWin = hanning(buffLen,'periodic') ;
% qDur = 4*buffLen/48000 ;

% mono audio files to be played
% fileName = 'E:\Work_shu\Project\CAMIL\the_CAMIL_dataset_part2\Ref\Emitted\0.wav';

% stero audio files to be played
% fileName = 'C:\Users\miya\Dropbox\computational auditory scene analysis\Drysound\sound_example_dry.wav';
fileName = '..\BinauralPlaybackTestSounds\NorahJonesSurround.wav';
% Surround audio files to be played
% fileName = 'F:\audioDatabase\SurroundRecordings\Avatar-04-0h29m42s-0h42m42s.wav';
% fileName = 'F:\audioDatabase\SurroundRecordings\VivaldiCheGiovaIlSospirar_24bit48kHz.wav';

% micrphone array recodering to be played
% fileName = 'F:\Dropbox\SphericalMircophone\20130926_Room1302\mic_092613_162130.wav';

% Reference distance for the "plane-waves"
% (the relative angles to the sources are calculated as if the sources sat
% on the surface of a sphere with radius refDist)
refDist = 2.0 ;

% Create file reader object and get number of channels
reader = dsp.AudioFileReader(fileName, ...
    'SamplesPerFrame',buffLen/2,'OutputDataType','single') ;
fileInfo = info(reader);
numChan = fileInfo.NumChannels;
% Surround channel mapping (remove sub channel if any)
if numChan == 1
    chanMap = 1;
elseif numChan == 2
    chanMap = 1:2;
elseif numChan == 5
    chanMap = 1:5 ;
elseif numChan == 6
    chanMap = [1 2 3 5 6] ;
elseif numChan == 32
    chanMap = 1:32;
else
    error('Wrong number of audio channels in the input file.')
end

% Load HRIRs
% load('HRTF\HRIR_ANL_ha') ;
% load(['C:\Users\miya\Documents\MATLAB\mcha\golay\results\HOAToolbox\HRTF\HRIR_ANL_' subjIni '.mat']) ;
load(['IRs\HRTF\HRIR_ANL_' subjIni '.mat']) ;
Hrir = hrir.impulseResponses ;
fsHrir = hrir.sampFreq ;
azHrir = hrir.sourceSphCoord(:,1) ;
elHrir = hrir.sourceSphCoord(:,2) ;
hrirLen = size(Hrir,1) ;
clear hrir

% Calculate HRTF magnitudes (or DTF)
HrtfMag = abs(fft(Hrir,2*buffLen)) ;
if dtf == true
    HrtfMag = exp(bsxfun(@minus,log(HrtfMag),mean(log(HrtfMag),3))) ;
end
HrtfMag = .99*HrtfMag/max(max(max(HrtfMag))) ;

% Minimum-phase filters
HrirMin = zeros(2*buffLen,2,length(azHrir)) ;
for I = 1 : 2
    for J = 1 : length(azHrir)
        [~,HrirMin(:,I,J)] = rceps(real(ifft(HrtfMag(:,I,J)))) ;
    end
end

% Calculate the delays between every original HRIR and its minimum-phase
% counterpart
numDir = length(azHrir) ;
numTap = size(Hrir,1) ;
HrirDel = zeros(numDir,2) ;
for I = 1 : numDir
    cor = xcorr(resample(Hrir(:,1,I),4,1),resample(HrirMin(:,1,I),4,1)) ;
    [~,maxIdx] = max(cor) ;
    HrirDel(I,1) = (maxIdx-8*buffLen)/4/fsHrir ;
    cor = xcorr(resample(Hrir(:,2,I),4,1),resample(HrirMin(:,2,I),4,1)) ;
    [maxCor,maxIdx] = max(cor) ;
    HrirDel(I,2) = (maxIdx-8*buffLen)/4/fsHrir ;
end

% Delay to be added to the interpolated Hrirs to avoid circular conv
% problems (if possible, depending on the buffer length)
cstDel = (buffLen-hrirLen)/4/48000 ;

% Reshape / normalise HRTF magnitudes
HrtfMag = .99*HrtfMag/max(max(max(HrtfMag))) ;
HrtfMag = permute(HrtfMag,[3 1 2]) ;

% Calculate spline parameters for HRTF magnitudes / delays 
HrtfMagSplParam = cell(2,1) ;
HrirDelSplParam = cell(2,1) ;
for I = 1 : 2
    HrtfMagSplParam{I} = ...
        SphericalSplineParam([azHrir,elHrir],HrtfMag(:,1:buffLen+1,I),1,1e-5) ;
    HrirDelSplParam{I} = ...
        SphericalSplineParam([azHrir,elHrir],HrirDel(:,I),1,1e-5) ;
end

% Surround speaker positions
if length(chanMap) == 1
    az0 = [0]' * pi/180 ;
    el0 = 0*pi/180*ones(size(az0)) ;
    [xyz0(:,1),xyz0(:,2),xyz0(:,3)] = sph2cart(az0,el0,refDist) ;
elseif length(chanMap) == 2
    az0 = [30 -30]' * pi/180 ;
    el0 = 0*pi/180*ones(size(az0)) ;
    [xyz0(:,1),xyz0(:,2),xyz0(:,3)] = sph2cart(az0,el0,refDist) ;
elseif length(chanMap) == 5
    % 5.1 surround sound layout
    az0 = [30 -30 0 110 -110]' * pi/180 ;
    el0 = 0*pi/180*ones(size(az0)) ;
    [xyz0(:,1),xyz0(:,2),xyz0(:,3)] = sph2cart(az0,el0,refDist) ;
elseif length(chanMap) == 32
    % "Plane-wave" basis
    numPlw = 32 ;
    [az0,el0] = SphericalTDesign(numPlw) ;
    [xyz0(:,1),xyz0(:,2),xyz0(:,3)] = sph2cart(az0,el0,refDist) ;
else
    error('Wrong number of audio channels in the input file.')
end


% HRTFs for the initial speaker positions
frq = (0:buffLen)'*48000/2/buffLen;
HrtfMagInt = zeros(length(chanMap),2*buffLen,2);
HrirDelInt = zeros(length(chanMap),2);
HrirMin = zeros(2*buffLen,length(chanMap),2);
HrtfInt = zeros(2*buffLen,length(chanMap),2);
for I = 1 : 2
    % HRTF magnitudes
    HrtfMagInt(:,1:buffLen+1,I) = ...
        SphericalSplineEval(HrtfMagSplParam{I},[az0,el0]) ;
    HrtfMagInt(:,buffLen+2:end,I) = ...
        HrtfMagInt(:,buffLen:-1:2,I) ;
    % Zero-Phase HRIRs
    HrirZer = real(ifft(HrtfMagInt,[],2)) ;
    % Hrir delays
    HrirDelInt(:,I) = cstDel + ...
        SphericalSplineEval(HrirDelSplParam{I},[az0,el0]) ;
    % Min-phase HRIRs
    for J = 1 : length(chanMap)
        [~,HrirMin(:,J,I)] = rceps(HrirZer(J,:,I).') ;
    end
    % Interpolated HRTFs
    HrtfInt(:,:,I) = fft(HrirMin(:,:,I)) ;
    HrtfInt(1:buffLen+1,:,I) = HrtfInt(1:buffLen+1,:,I) ...
        .* exp(-1i*2*pi*frq*HrirDelInt(:,I).') ;
    HrtfInt(buffLen+2:end,:,I) = ...
        conj(HrtfInt(buffLen:-1:2,:,I)) ;
end

% Load diffuse reverb filters,
load stereoRirListeningRoom.mat;

% Convolve the RIRs with diffuse HRIRs -> BRIRs
% if dtf == false
%     load(['E:\Work_shu\Project\matlabG4interface\matlabScripts\diffuseHRIRs\diffuseHRIRs_' subjIni '.mat']) ;
%     rir = fftfilt(diffHrir,rir) ;
% end

% divide BRIRs in chunks and calculate the FFT for each chunk
numChunks = ceil(size(rir,1)/buffLen) ;
rir = [rir ; zeros(numChunks*buffLen-size(rir,1),2)] ;
rirFft = zeros(2*buffLen,numChunks,2) ;
for I = 1 : numChunks
    idx = (I-1)*buffLen + (1:buffLen) ;
    for J = 1 : 2
        rirFft(:,I,J) = fft(rir(idx,J),2*buffLen) ;
    end
end

% Normalise energy of the reverb to get an appropriate DRR
drr = 15 ;
energDir = sum(sum(sum(real(ifft(HrtfInt)).^2))) ;
energRvb = sum(sum(rir.^2)) ;
rirFft = rirFft * 10^(-drr/20) * sqrt(energDir/energRvb) ;

% Create audio player object
player = dsp.AudioPlayer('BufferSizeSource','Property','BufferSize',buffLen);

% Create audio file writer
filewriter = dsp.AudioFileWriter('Binaural_Static.wav','SampleRate',reader.SampleRate,'FileFormat','WAV');

% Load Head tracker data
% load 'RecPosOri'
load 'SimPosOri_Pan'

% initialise position
L = length(data);
posIni = data(1,1:3);
timCur = data(1,7);
i = 0;

%% Filter signals and play back
binBuff = zeros(2*buffLen,2);
inpBuff = zeros(buffLen,numChan);
spkBuff = zeros(buffLen,length(chanMap)) ;
omniFftBuff = zeros(2*buffLen,numChunks) ;
tic
while ~isDone(reader) % && i < L
    % New chunk of 5.1 signals
    inpBuff(end/2+1:end,:) = step(reader) ;
    % FFT of the new chunk of speaker signals
    spkBuff = bsxfun(@times,buffWin,10^(10/20)*inpBuff(:,chanMap)) ;
    spkFft = fft(spkBuff,2*buffLen) ;
    % Get tracker position/orientation and update rotation matrix
    if HT == 1
        pos = data(i+1,1:3);
        pos = pos/100 ; % convert to meter
        eulerAng = data(i+1,4:6);
        timeStep = data(i+1,7)-timCur;
        i = mod(i + 1, L);
        if toc > Tim
            tic
            % Calculate rotation matrix
            Rx = [ 1                0                 0 ;
                   0 cos(eulerAng(3)) -sin(eulerAng(3)) ;
                   0 sin(eulerAng(3))  cos(eulerAng(3)) ];
            Ry = [ cos(eulerAng(2)) 0 sin(eulerAng(2)) ;
                    0               1                0 ;
                   -sin(eulerAng(2)) 0 cos(eulerAng(2)) ];
            Rz = [ cos(eulerAng(1)) -sin(eulerAng(1)) 0 ;
                    sin(eulerAng(1))  cos(eulerAng(1)) 0 ;
                    0                 0                1 ];
            R = Rz*Ry*Rx ;
            % Update speaker positions
            xyz = bsxfun(@minus,xyz0,pos) ;
            xyz = xyz * R ;
            [azSpk,elSpk] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3)) ;
            % Interpolate HRTFs
            for I = 1 : 2
                % HRTF magnitudes
                HrtfMagInt(:,1:buffLen+1,I) = ...
                    SphericalSplineEval(HrtfMagSplParam{I},[azSpk,elSpk]) ;
                HrtfMagInt(:,buffLen+2:end,I) = ...
                    HrtfMagInt(:,buffLen:-1:2,I) ;
                % Hrir delays
                HrirDelInt(:,I) = cstDel + ...
                    SphericalSplineEval(HrirDelSplParam{I},[azSpk,elSpk]) ;
                % Zero-Phase HRIRs
                HrirZer = real(ifft(HrtfMagInt,[],2)) ;
                % Min-phase HRIRs
                for J = 1 : length(chanMap)
                    [~,HrirMin(:,J,I)] = rceps(HrirZer(J,:,I).') ;
                end
                % Interpolated HRTFs
                HrtfInt(:,:,I) = fft(HrirMin(:,:,I)) ;
                HrtfInt(1:buffLen+1,:,I) = HrtfInt(1:buffLen+1,:,I) ...
                    .* exp(-1i*2*pi*frq*HrirDelInt(:,I).') ;
                HrtfInt(buffLen+2:end,:,I) = ...
                    conj(HrtfInt(buffLen:-1:2,:,I)) ;
            end
        end
    end
    
    % FFT of the binaural signal buffer
    binFft = squeeze(sum(bsxfun(@times,HrtfInt,spkFft),2)) ;
    
    % Add reverb
    omniFftBuff(:,1) = mean(spkFft,2) ;
    binFft = binFft + squeeze(sum(bsxfun(@times,rirFft,omniFftBuff),2)) ;
  
    % Add new binaural signals to the end of the previous buffer
    binBuff = binBuff + real(ifft(binFft)) ;
    
    % Play back binaural signals
    step(player,binBuff(1:buffLen/2,:));    
    step(filewriter, binBuff(1:buffLen/2,:));
    
    % Shift buffers
    inpBuff = [inpBuff(end/2+1:end,:);zeros(buffLen/2,numChan)] ;
    binBuff = [binBuff(buffLen/2+1:end,:);zeros(buffLen/2,2)] ;
    omniFftBuff = [zeros(2*buffLen,1),omniFftBuff(:,1:end-1)] ;
end

% Release player/reader
release(reader)
release(player)
% release(filewriter)
%%
