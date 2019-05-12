close all, clear sound
%Bruno E. Gracia Villalobos
%EE 3523 - Discrete Signals and Systems
%Final Project (5/10/19)
%Professor Grigoryan

disp("> Bruno E. Gracia Villalobos - Signal Convolution with Filters");
disp("> 5/3/19");
disp("> EE 3523");
disp("> Professor Grigoryan");

%~~~~~~~~~~~~~~~CONFIG PARAMETERS~~~~~~~~~~~~~~~~~~~%

duration = 40;          %how long should signal be (seconds)
begin = 0;              %reference point to start (seconds)
chunkSize = 10;     %how many seconds per chunk to do convolution
filePath = "mike.wav"; %file path to .wav relative to home directory

%~~~~~~~~~~~~~~~CONFIG PARAMETERS~~~~~~~~~~~~~~~~~~~%

%Read in info from WAV file
[samples, fs] = audioread(filePath);
info = audioinfo(filePath);
disp(info);

%Get info from file
channels = info.NumChannels;
sampleRate = info.SampleRate;
timeLimit = info.Duration;

%specify beginning sample @ 1 due to  MATLAB constraints
if (sampleRate*begin) == 0
    lowerBound = 1;
else
    lowerBound = ceil(sampleRate*begin);
end

%check if user specified bounds are possible given wav size
if(begin+duration) >timeLimit
    upperBound = ceil(timeLimit*sampleRate);
else
     upperBound = lowerBound + ceil(duration*sampleRate);
end

%specify how much of the .wav file to use
if channels == 2
    lchannel = samples( lowerBound : upperBound-1 ,1) .' ;
    rchannel = samples( lowerBound : upperBound-1 ,2) .' ;
else
   lchannel = samples( lowerBound : upperBound-1).' ;
   rchannel = lchannel;
end
 
%from now on, refer to total samples to length of our wav sample
numSamples = length(lchannel);

%Assemble both channels into a track for playback
songSnippet = [ lchannel.' rchannel.'];

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

%~~~~~~~~~~~~~~~~~~~~~~~~PART 2A~~~~~~~~~~~~~~~~~~~~~~
%                                       x(n)*h(n) Convolution by sum
tic

fprintf("%s\n> ", "Direct Convolution");
%LEFT channel
convLeft = lchannel;
convFilter = filter(length(filter) : -1 : 1); %mirror

%specify bounds of convolution for for loop
filterLow = ceil(length(filter)/2);
filterHigh = floor(length(filter)/2);

for n=filterLow:numSamples-filterHigh
    p=lchannel(n-filterHigh:n+filterHigh);
    convLeft(n)=sum(p.*convFilter);
end

%RIGHT channel
convRight = rchannel;

for n=filterLow:numSamples-filterHigh
    p=rchannel(n-filterHigh:n+filterHigh);
    convRight(n)=sum(p.*convFilter);
end

toc
fprintf("\n");

%Plot Filter
figure;
plot(filterAxis, filter, 'r');
title("Filter plot");

printName("");

% PLOT L AND R CHANNELS
%Create timeAxis from begin time to end time.
%Subtract 1/sampleRate because MATLAB adds 1
timeAxis = 0 : (1/sampleRate) : ( length(lchannel)/sampleRate - 1/sampleRate );

%offset timeAxis depending on begin time
timeAxis = timeAxis + begin;
%tAxis = begin:(1/sampleRate): (begin+ (chunkSamples/sampleRate) + ((length(filter)-2)/sampleRate));

%Plot the signal sample from the .wav without filtering
if channels == 2
    figure;
    subplot(2,2,1)
    plot(timeAxis,lchannel, 'b');
    title("Left Channel");

    subplot(2,2,2)
    plot(timeAxis,rchannel, 'g');
    title("Right Channel");

    %Plot LChannel Convolved with Filter
    subplot(2,2,3);
    plot(timeAxis,convLeft, 'b');
    title("L CHANNEL convolved with FOR LOOP");

    %Plot RChannel Convolved with Filter
    subplot(2,2,4);
    plot(timeAxis,convRight, 'g');
    title("R CHANNEL convolved with FOR LOOP");
else
    figure;
    subplot(2,1,1)
    plot(timeAxis,lchannel, 'b');
    title("Signal");
    
    subplot(2,1,2)
    plot(timeAxis,convLeft, 'g');
    title("Signal with Direct Convolution");
end

printName("");

%PLAY LPF song
songLPF = [convLeft.' convRight.'];
%sound(songLPF, sampleRate);

%~~~~~~~~~~~~~~~~~~~~~~~~CONVOLUTION BY FFT~~~~~~~~~~~~~~~~~~~~~
tic

fprintf("%s\n> ", "FFT Convolution");

chunkSamples = sampleRate*chunkSize; %set the size of each chunk
chunksTotal = ceil(duration/chunkSize); %calculate how many chunks are needed
totChunkSamples = chunksTotal*chunkSamples; %calculate how many samples total for all chunks

%find out the bounds of each chunk in samples
chunks = (begin*sampleRate) : chunkSamples : (begin*sampleRate)+totChunkSamples;
remSamples = totChunkSamples-length(lchannel); %zero padding

%if there is zero padding needed, append channels with zeros
if remSamples>0
    lchannel = [lchannel zeros(1, remSamples)];
    rchannel = [rchannel zeros(1, remSamples)];
end

chunks = chunks + 1;    %to adhere to MATLAB guidelines
chunksTotal = length(chunks) - 1; %constant for total chunks

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
fprintf("\n");

%resize timeAxis to account for the added samples from convolution
timeAxis = 0:(1/sampleRate):( length(conv_l) * (1/sampleRate) - 1/sampleRate);
%add again offset to display plot
timeAxis = timeAxis + begin;

%PLOT RESULTS

if channels ==2
    figure;
    subplot(2,1,1);
    plot(timeAxis, conv_l, 'b');
    title("L CHANNEL Convolved with FFT");

    subplot(2,1,2);
    plot(timeAxis, conv_r, 'g');
    title("R CHANNEL Convolved with FFT");

else
    figure;
    subplot(2,1,1)
    plot(timeAxis,[lchannel zeros(1, length(filter)-1)], 'b');
    title("Signal");
    
    subplot(2,1,2)
    plot(timeAxis,conv_l, 'g');
    title("Signal with FFT Convolution");
end
printName("");


%sound([conv_l.' conv_r.'], sampleRate);

%~~~~~~~~~~~~~~~~~~~~~~~~PART 2B-1~~~~~~~~~~~~~~~~~~~~~~
tic

fprintf("%s\n> ", "Overlap Add Convolution - Convolving each Chunk");
%set aside memory for overlap add method.
chunkV_l = zeros(chunksTotal, chunkSamples + (length(filter)-1) );
chunkV_r = chunkV_l;

%remove begin offset to use as indexes for convolution:
if (chunks(1)-lowerBound) < 1
    chunks = chunks - lowerBound+1;
else
    chunks = chunks - lowerBound;
end

%make loop to make a matrix to hold the chunks
for i=1:(chunksTotal)
    %get the chunks from song for each channel
    %assign to the first ten seconds the fetched samples
    chunkV_l(i,1:chunkSamples) = lchannel(chunks(i):chunks(i+1)-1);;
    chunkV_r(i,1:chunkSamples) = rchannel(chunks(i):chunks(i+1)-1);
end

%%%%%%%%%%setup variables for FFT%%%%%%%%%%%%%

%same as chunkV_r setup but for storing FFTs
chunkV_r_fft = zeros(chunksTotal, length(chunkV_l(1,:)));
chunkV_l_fft = chunkV_r_fft;

%get the frequency bins and store in another matrix
for i=1:chunksTotal
    chunkV_l_fft(i,:) = fft(chunkV_l(i,:));
    chunkV_r_fft(i,:) = fft(chunkV_r(i,:));
end

%now resize the filter vector for each tensec chunk
fft_f = fft(filter,length(chunkV_l_fft(1,:)));

%now that both the filter and  chunk are the same size
%with zero padding, time to convolve!

%setup again matrices to store the convolutions
overlap_fft_l = zeros(chunksTotal, length(chunkV_l(1,:)));
overlap_fft_r = overlap_fft_l;

% F(signal)*F(filter) for each channel
for i=1:chunksTotal
    overlap_fft_l(i,:) = chunkV_l_fft(i,:) .* fft_f; 
    overlap_fft_r(i,:) = chunkV_r_fft(i,:) .* fft_f; 
end

%now time to go back to time domain to play the filtered song!
overlap_l = zeros(chunksTotal, length(chunkV_l(1,:)));
overlap_r = overlap_l;

%inverse fft for each chunk and store in l and r matrices
for i=1:chunksTotal
    overlap_l(i,:) = ifft( overlap_fft_l(i,:) );
    overlap_r(i,:) = ifft( overlap_fft_r(i,:) );
end

toc
fprintf("\n");

%~~~~~~~~~~~~~~~~~~~~~~~~LOOP TO PRINT FIGURES~~~~~~~~~~~~~~~~~~~~~~
%create a time axis to represent each chunk in time
%(chunkSamples/sampleRate) = length of chunk in time
%(length(filter)-2)/sampleRate)) = length of filter in time
tAxis = begin:(1/sampleRate): (begin+ (chunkSamples/sampleRate) + ((length(filter)-2)/sampleRate));

%loop counter
currentChunk =1;

%add again begin offset for plotting
chunks = chunks + lowerBound-1;

if(channels == 1)
   for i=1:chunksTotal
       
       %if 4 plots have been printed, create new figure to keep organized
       if (i==1) || (mod(i,4) ==1 )
            figure;
            printName("");
       end
       
       %caption for plot
       caption =   "[" + chunks(currentChunk) + ", " + (chunks(currentChunk+1)-1) + "] samples";
        
       plotNum = mod(i,4);
       if plotNum == 0
           plotNum = 4;
       end
       %plot channel figure
       subplot(4,channels ,plotNum);
       plot(tAxis, overlap_l(currentChunk, :), 'b');
       title(caption);
       
       %offet for next loop run
       tAxis = tAxis + chunkSize;
       currentChunk = currentChunk + 1;
        
   end
else
    for i=1 : 2 : (chunksTotal*channels)
        %generate a new figure if 4 rows have been displayed
        if (i==1) || (mod(i,8)==1)
            %fprintf("new figure");
            figure;
            printName("");
        end

        %diagnostics
        counter = "i is " +i;
        ctrChunk = "currentChunk is " + currentChunk;
        
        %disp(counter);
        %disp(ctrChunk);
        
        %figure labels
        leftCaption =   "[" + chunks(currentChunk) + ", " + (chunks(currentChunk+1)-1) + "] samples of L channel";
        rightCaption = "[" + chunks(currentChunk) + ", " + (chunks(currentChunk+1)-1) + "] samples of R channel";
        
        plotNum = mod(i,8);
        if plotNum == 0
            plotNum = 8;
        end
       
        %plot left channel
        subplot(4,channels, plotNum);
        plot(tAxis, overlap_l(currentChunk, :), 'b');
        title(leftCaption);

        %plot right channel
        subplot(4,channels , plotNum+1);
        plot(tAxis, overlap_r(currentChunk, :), 'g');
        title(rightCaption);

        %offet for next loop run
        tAxis = tAxis + chunkSize;
        currentChunk = currentChunk + 1;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~PART 2B-2~~~~~~~~~~~~~~~~~~~~~~
tic

fprintf("%s\n> ", "Overlap Add Convolution - Overlap and Add");
%setup final time domain vectors for each channel
fft_conv_l = zeros(1, length(lchannel) + length(filter)-1);
fft_conv_r = fft_conv_l;

overlapLen = length(overlap_fft_l(1,:));

%now add and overlap
for i=1:chunksTotal-1
    if i==1
        fft_conv_l(i:overlapLen) = overlap_l(i,:);
        fft_conv_r(i:overlapLen) = overlap_r(i,:);
    end
    
    %here, each 10 sec chunk is being added and between 480000 and 480020,
    %overlap is being added.
    fft_conv_l(i*chunkSamples:i*chunkSamples+overlapLen-1) = ...
        fft_conv_l(i*chunkSamples:i*chunkSamples+overlapLen-1) ...
        + overlap_l(i+1,:);
    
    %do the same for the right channel
    fft_conv_r(i*chunkSamples:i*chunkSamples+overlapLen-1) = ...
        fft_conv_r(i*chunkSamples:i*chunkSamples+overlapLen-1) ...
        + overlap_r(i+1,:);  
end

toc
fprintf("\n");

%resize timeAxis to fit overlap add convolution length
timeAxis = 0:(1/sampleRate):( length(fft_conv_l) * (1/sampleRate) - 1/sampleRate);

%offset again for begin time
timeAxis = timeAxis + begin;

if channels == 2
    figure;
    subplot(2,1,1);
    plot(timeAxis, fft_conv_l, 'b');
    title("L channel signal after OVERLAP ADD");

    subplot(2,1,2);
    plot(timeAxis, fft_conv_r,  'g');
    title("R channel signal after OVERLAP ADD");

else
    figure;
    subplot(2,1,1)
    plot(timeAxis,[lchannel zeros(1, length(filter)-1)], 'b');
    title("Signal");
    
    subplot(2,1,2)
    plot(timeAxis,fft_conv_l, 'g');
    title("Signal with Overlap Add Convolution");
end
%play filtered song
G = [fft_conv_l.' fft_conv_r.'];
sound( G , sampleRate);

%calculate MSR error
MSR_L = calcMSR(conv_l, fft_conv_l);
MSR_R = calcMSR(conv_r, fft_conv_r);

if channels == 2
    fprintf('MSR error for Left Channel = %12.10f \n',MSR_L);
    fprintf('MSR error for Right Channel = %12.10f \n',MSR_R);
else
    fprintf('MSR error for Signal = %12.10f \n',MSR_L);
end


printName("");

%~~~~~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~%

function printName(~)
    %legend('B.E.G.V.','location','EastOutside')
    MyBox = uicontrol('style','text');
    set(MyBox,'String','B.E.G.V.');
    set(MyBox,'Position',[20,20,70,20]);
    set(MyBox,'BackgroundColor','green');
end

function MSR = calcMSR(sig1, sig2)
    MSR = sqrt(sum((sig1 - sig2).^2 )) / length(sig1);
end