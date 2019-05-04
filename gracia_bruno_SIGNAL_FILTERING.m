close all
%Bruno E. Gracia Villalobos
%EE 3523 - Discrete Signals and Systems
%Final Project
%Professor Grigoryan

%Read in my song "ascend - bruno G*"
% "native" => samples are 24 bits stored as int32's
[samples, fs] = audioread("ascend_1604.wav", "double");
info = audioinfo("ascend_1604.wav");

%Read in info from WAV file
disp("> Bruno E. Gracia Villalobos - Signal Convolution with Filters");
disp("> 5/3/19");
disp("> EE 3523");
disp("> Professor Grigoryan");
disp(info);
%{

duration             = info.Duration; %in seconds
numSamples       = info.TotalSamples;

disp(duration);
%}

%Get L and R Channels from WAV and make into row
%get 40 seconds of my song
sampleRate = info.SampleRate;
begin = 28; %start at second 28
duration = 40;

lchannel = samples( (sampleRate*begin) : ...
    (sampleRate*(begin+duration))-1 ,1) .' ;

rchannel = samples( (48000*begin) : ...
    (sampleRate*(begin+duration))-1 ,2) .' ;

numSamples = length(lchannel);

%Assemble both channels into a track for playback
songSnippet = [ lchannel.' rchannel.'];

%play the audio tracks
sound(songSnippet,sampleRate);

%Create time axis for each individual sample
%timeVector  = linspace(0, duration, numSamples);
temp = begin:(1/sampleRate):(begin+duration);
timeAxis = temp(1:numSamples);

%Create Triangle Filter filter(11) is the center
filter = [1 2 3 3 4 4 5 6 6 7 9 7 7 6 6 5 5 4 4 3 1];
%{
filter = [-659 -1915 -2005 -358 1679 1089 -1853 -2807 ...
    2077 10186 14235 10186 2077 -2807 -1853 1089 ...
1679 -358 -2005 -1915 -659];
%}
%filter = [3 -1 -3 -1 -7 2 25 -3 -44 1 53 1 -44 -3 25 2 -7 -1 -3 -1 3];
%BPF 1-2KHZ 
%filter = [0 0 2 0 -6 0 10 0 -15 0 16 0 -15 0 10 0 -6 0 2 0 0];

k = sum(filter); % normalize
filter = filter/k;
filterAxis = 1:1:length(filter);

%%%%%%%%%%%%%%%%%PART A%%%%%%%%%%%%%
%%%%%% x(n)*h(n) Convolution by sum
tic
%LEFT channel
convLeft = lchannel;
convFilter = filter(length(filter) : -1 : 1); %mirror

for n=11:numSamples-10
    p=lchannel(n-10:n+10);
    convLeft(n)=sum(p.*convFilter);
end

%RIGHT channel
convRight = rchannel;

for n=11:numSamples-10
    p=rchannel(n-10:n+10);
    convRight(n)=sum(p.*convFilter);
end
toc

%Plot Filter
figure;
plot(filterAxis, filter, 'r');
title("Filter of 21 taps");

MyBox = uicontrol('style','text');
set(MyBox,'String','B.E.G.V.');
set(MyBox,'Position',[20,20,70,20]);
set(MyBox,'BackgroundColor','green');

% PLOT L AND R CHANNELS
figure;
subplot(2,2,1)
plot(timeAxis,lchannel, 'b');
title("Left Channel (from 00:28:00 to 00:68:00)");

subplot(2,2,2)
plot(timeAxis,rchannel, 'g');
title("Right Channel (from 00:28:00 to 00:68:00)");

%Plot LChannel Convolved with Filter
subplot(2,2,3);
plot(timeAxis,convLeft, 'b');
title("L CHANNEL convolved with FOR LOOP");

%Plot RChannel Convolved with Filter
subplot(2,2,4);
plot(timeAxis,convRight, 'g');
title("R CHANNEL convolved with FOR LOOP");

MyBox = uicontrol('style','text');
set(MyBox,'String','B.E.G.V.');
set(MyBox,'Position',[20,20,70,20]);
set(MyBox,'BackgroundColor','green');

%PLAY LPF song
songLPF = [convLeft.' convRight.'];
%sound(songLPF, sampleRate);

%%%%%%%%%%%%%%%%%Convolution by FFT%%%%%%%%%%%%%
tic
%Setup FFT of L and R Channels, and Filter
fft_l = fft( [ lchannel zeros(1, length(filter)-1) ] );
fft_r = fft( [ rchannel zeros(1, length(filter)-1) ] );
%fft_f = fft([filter zeros(1, length(lchannel) - length(filter)) ]); %not
%needed, automatic padding
fft_f = fft(filter, length(fft_l));

%Convolve in frequency domain
convF_l = fft_l .* fft_f;
convF_r = fft_r .*fft_f;

%Go back to time domain
conv_l = ifft(convF_l);
conv_r = ifft(convF_r);

toc
%resize timeAxis to account for the added samples from convolution
timeAxis =  begin:(1/sampleRate):(begin+duration + (19/sampleRate));

figure;
subplot(2,1,1);
plot(timeAxis, conv_l, 'b');
title("L CHANNEL Convolved with FFT");

subplot(2,1,2);
plot(timeAxis, conv_r, 'g');
title("R CHANNEL Convolved with FFT");

%sound([conv_l.' conv_r.'], sampleRate);

%%%%%%%%%%%%PART B-1%%%%%%%%%%%%%%%%%%%
tic
tenSec = sampleRate*10; %constant
chunks = 0: (tenSec): numSamples; %mark the bounds of each chunk
chunks = chunks + 1;
totChunks = length(chunks) - 1; %constant for total chunks

%set aside memory for overlap add method.
%this is a 4 row by 480000 + 20 zeros vector for each convolution
chunkV_l = zeros(totChunks, tenSec + (length(filter)-1) );
chunkV_r = chunkV_l;

%make loop to make a matrix to hold 10 second chunks
for i=1:totChunks
    %get the 10 sec chunks from song for each channel
    L = lchannel(chunks(i):chunks(i+1)-1);
    R = rchannel(chunks(i):chunks(i+1)-1);
    
    %assign to the first ten seconds the fetched samples
    chunkV_l(i,1:tenSec) = L;
    chunkV_r(i,1:tenSec) = R;
end

%%%%%%%%%%setup variables for FFT%%%%%%%%%%%%%

%same as chunkV_r setup but for storing FFTs
chunkV_r_fft = zeros(totChunks, length(chunkV_l(1,:)));
chunkV_l_fft = chunkV_r_fft;

%get the frequency bins and store in another matrix
for i=1:totChunks
    chunkV_l_fft(i,:) = fft(chunkV_l(i,:));
    chunkV_r_fft(i,:) = fft(chunkV_r(i,:));
end

%now resize the filter vector for each tensec chunk
fft_f = fft(filter,length(chunkV_l_fft(1,:)));

%now that both the filter and 10 sec chunk are the same size
%with zero padding, time to convolve!

%setup again matrices to store the convolutions
overlap_fft_l = zeros(totChunks, length(chunkV_l(1,:)));
overlap_fft_r = overlap_fft_l;

% F(signal)*F(filter) for each channel
for i=1:totChunks
    overlap_fft_l(i,:) = chunkV_l_fft(i,:) .* fft_f; 
    overlap_fft_r(i,:) = chunkV_r_fft(i,:) .* fft_f; 
end

%now time to go back to time domain to play the filtered song!
overlap_l = zeros(totChunks, length(chunkV_l(1,:)));
overlap_r = overlap_l;

%inverse fft for each chunk and store in l and r matrices
for i=1:totChunks
    overlap_l(i,:) = ifft( overlap_fft_l(i,:) );
    overlap_r(i,:) = ifft( overlap_fft_r(i,:) );
end

toc
%PLOT
temp_t = begin:(1/sampleRate):(begin+10 + (19/sampleRate));
tAxis = temp_t(1:length(temp_t));

figure;
subplot(4,2,1);
plot(tAxis, overlap_l(1,:), 'b');
title("1 to 480020th samples of L channel");

subplot(4,2,2);
plot(tAxis, overlap_r(1,:), 'g');
title("1 to 480020th samples of R channel");

temp_t = (begin+10):(1/sampleRate):(begin+20 + (19/sampleRate));
tAxis = temp_t(1:length(temp_t));

subplot(4,2,3);
plot(tAxis, overlap_l(2,:), 'b');
title("480020th to 980040th samples of L channel");

subplot(4,2,4);
plot(tAxis, overlap_r(2,:), 'g');
title("480020th to 980040th  samples of R channel");

temp_t = (begin+20):(1/sampleRate):(begin+30 + (19/sampleRate));
tAxis = temp_t(1:length(temp_t));

subplot(4,2,5);
plot(tAxis, overlap_l(3,:), 'b');
title("980040th to 1440060th samples of L channel");

subplot(4,2,6);
plot(tAxis, overlap_r(3,:), 'g');
title("980040th to 1440060th  samples of R channel");

temp_t = (begin+30):(1/sampleRate):(begin+40 + (19/sampleRate));
tAxis = temp_t(1:length(temp_t));

subplot(4,2,7);
plot(tAxis, overlap_l(4,:), 'b');
title("1440060th to 1920080th samples of L channel");

subplot(4,2,8);
plot(tAxis, overlap_r(4,:), 'g');
title("1440060th to 1920080th  samples of R channel");


%%%%%%%%%%%%PART B-2%%%%%%%%%%%%%%%%%%%

tic
%setup final time domain vectors for each channel
fft_conv_l = zeros(1, length(lchannel) + length(filter)-1);
fft_conv_r = fft_conv_l;

overlapLen = length(overlap_fft_l(1,:));


%now add and overlap
for i=1:totChunks-1
    if i==1
        fft_conv_l(i:overlapLen) = overlap_l(i,:);
        fft_conv_r(i:overlapLen) = overlap_r(i,:);
    end
    
    %here, each 10 sec chunk is being added and between 480000 and 480020,
    %overlap is being added.
    fft_conv_l(i*tenSec:i*tenSec+overlapLen-1) = ...
        fft_conv_l(i*tenSec:i*tenSec+overlapLen-1) ...
        + overlap_l(i+1,:);
    
    %do the same for the right channel
    fft_conv_r(i*tenSec:i*tenSec+overlapLen-1) = ...
        fft_conv_r(i*tenSec:i*tenSec+overlapLen-1) ...
        + overlap_r(i+1,:);  
end

toc

figure;
subplot(2,1,1);
plot(timeAxis, fft_conv_l, 'b');
title("L channel signal after OVERLAP ADD");

subplot(2,1,2);
plot(timeAxis, fft_conv_r,  'g');
title("R channel signal after OVERLAP ADD");


%play filtered song
G = [fft_conv_l.' fft_conv_r.'];
%sound( G , sampleRate);

%calculate MSR error
MSR_L = sqrt(sum((conv_l - fft_conv_l).^2 )) / length(conv_l);
MSR_R = sqrt(sum((conv_r - fft_conv_r).^2 )) / length(conv_r);

fprintf('  MSR error for Left Channel = %6.4f \n',MSR_L);
fprintf('  MSR error for Right Channel = %6.4f \n',MSR_R);


%legend('B.E.G.V.','location','EastOutside')

MyBox = uicontrol('style','text');
set(MyBox,'String','B.E.G.V.');
set(MyBox,'Position',[20,20,70,20]);
set(MyBox,'BackgroundColor','green');

%CLEANUP
%clear sound

