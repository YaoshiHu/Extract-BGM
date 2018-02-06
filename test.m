%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Signal Processing Final Project- Columbia University           %
%            Extract Background Music from Soundtracks                   %
%           Minglei Gu, Leying Hu, Yaoshi Hu, Xingzhi Li                 %
%          {mg3847, lh2871, yh2950, xl2680}@columbia.edu                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% The main function for signal processing
function y = test(Data,fs)

    % Collect sizing from the input
    [SampleSize, ChannelSize] = size(Data);

    % Time length for every sample
    SampleSizeInSeconds = fs/SampleSize; 

    % Find the min num larger than sampling time 
    % that is the exponential of 2 for FFT
    N = 2.^nextpow2(fs*SampleSizeInSeconds);
    
    % Generate STFT for both channels
    X = [];
    waitbar1 = waitbar(0,'Start Converting Signals...');
    for i = 1:ChannelSize
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stft Function Reference                                       %
        % Matlab.com. (2017). Inverse Short-Time Fourier Transformation %
        % (ISTFT) with Matlab Implementation - File Exchange -          %
        % MATLAB Central. [online] Mathworks.com.                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xi = stft(Data(:,i), N, N/4, N*2, fs);
        % Compile FT results of the channels separately
        % X dimensions: Freq Response, Time Variance, Channel Num
        X = cat(3,X,Xi);
        waitbar(i/ChannelSize,waitbar1)
    end
    close(waitbar1)
    % Extract left half of the symmetric signal
    V = abs(X(1:N/2+1,:,:));
    
    % Find periods of the sinal
    p = correlation3D(V);

    %Ensure the output audio file length is not changing
    y = zeros(SampleSize,ChannelSize);                                
    waitbar2 = waitbar(0,'Outputting Signals...');
    for i = 1:ChannelSize
        % Finding similarities from periodic intervals
        Mi = similarsample(V(:,:,i),p);
        % Mirror back the signal
        Mi = cat(1,Mi,flipud(Mi(2:end,:)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % istft Function Reference                                      %
        % Matlab.com. (2017). Inverse Short-Time Fourier Transformation %
        % (ISTFT) with Matlab Implementation - File Exchange -          %
        % MATLAB Central. [online] Mathworks.com.                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        yi = istft(Mi.*X(:,:,i), N, N/4, N*2, fs);

        if length(yi) < SampleSize
            % Zero padding the length
            yi = [yi, zeros(1, SampleSize-length(yi))];
        end
        y(:,i) = yi(1:SampleSize); 
        waitbar(i/ChannelSize,waitbar2)
    end
    close(waitbar2)
    % Filter out the noise of the output signal
    y = ProcessFIR(y,10000,fs);
end

%%
% Finding similarities from periodic intervals
function simi = similarsample(signal,period)
    [m,n] = size(signal);
    % Number of repeating BGM in the signal
    repeat = ceil(n/period); 
    % Number of samples for the remaining signal is shorter than 1 period
    % so the padding zeros is needed
    lastlength = n-(repeat-1)*period; 
    template = zeros(m,period);
    for i = 1:period
        if(i<=lastlength)
            for j = 1:repeat
                template(:,i) = template(:,i)+signal(:,i+period*(j-1));
            end
            template(:,i) = template(:,i)/repeat;
        else
            for j = 1:repeat-1
                template(:,i) = template(:,i)+signal(:,i+period*(j-1));
            end
            template(:,i) = template(:,i)/(repeat-1);
        end
    end
    % Forge a signal composed of duplicated medium repeating sample
    forgeSignal = repmat(template, [1,repeat]); 
    forgeSignal = forgeSignal(:,1:n);
    % Magnitude cannot exceed original
    forgeSignal = min(forgeSignal, signal); 
    % In case the original value equal to 0 and output NAN
    simi = forgeSignal./(signal+eps); 
end

%%
% Input should have at least 200 points or this will output a wrong answer
function index = correlation3D(input)
    a = mean(input.^2,3);
    [m,n] = size(a);
    similarity = zeros(m,n);
    waitbar1 = waitbar(0,'Computing Correlation...');
    for i = 1:m
        similarity(i,:) = correlation(a(i,:));
        waitbar(i/m,waitbar1)
    end
    close(waitbar1)
    y = zeros(1,n);
    for i = 1:n
        y(i) = sum(similarity(:,i));
    end
    y = y/m;
    y(1) = [];
    index = 0;
    % Only count from the 2nd period sample
    % since the first period has the largest similarity with itself
    for i = ceil(n/50):length(y)-1 
         % Find the local maximum of all the similarity and
         % find the period of the song
        if(y(i)>y(i-1)&&y(i)>y(i+1)&&y(i)>0.8*max(y(ceil(n/50):end)))
            index = i;
            break
        end
    end
end

function y = correlation(a)
    % Set different weight of different frequency according to their value
    % Normally, the first few seconds will be pure background music
    b = a(1:ceil(length(a)/50));  
    c = xcorr2(a,b);
    bmod = sqrt(sum(b(:).^2));
    amod = zeros(1,length(a)+length(b)-1);
    % Adding zeros before and after sequence a to compute the normalized value
    tmp = [zeros(1,length(b)-1),a,zeros(1,length(b)-1)];  
    for i = 1:length(amod)
        amod(i) = sqrt(sum(tmp(i:i+length(b)-1).^2));
    end
    y = c./amod/bmod;
    y(1:length(b)-1) = [];
end

%%
% Filter out the noise of the output signal
function Background = ProcessFIR(Data, Wc, fs)
    
    % Process Parks-McClellan optimal FIR filter order estimation
   
    % Frequency band edges ranges [0, fs/2]
    f = [Wc, Wc*2];
    % Desired amp for the output (input signals have amp from 0 to 1)
    a = [1 0];
    % Max allowable deviation or ripples (Trials and Errors result)
    dev = [0.01 0.001];
    % n->approximate order
    % f0->normalized freq band edges
    % a0->freq band amp
    % w->weights that meet input specifications
    [n,f0,a0,w] = firpmord(f/(fs/2), a, dev);
    
    %Process Parks-McClellan optimal FIR low press filter design
    %order n FIR filter coefficients
    y = firpm(n, f0, a0, w);
    Background = filter(y, 1, Data);
end