clc
clear
close all

%%
% Input by the User
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in the name of the audio file
input_name = 'wonderful_world_input.wav';

% Fill in how many times to repeat the correlation process
N = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% run_me.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital Signal Processing Final Project- Columbia University           %
%            Extract Background Music from Soundtracks                   %
%           Minglei Gu, Leying Hu, Yaoshi Hu, Xingzhi Li                 %
%          {mg3847, lh2871, yh2950, xl2680}@columbia.edu                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Define input and output file locations
wav_file = input_name;
background_output = ['background_', input_name(1:strfind(input_name, '.') -1), '_processed.wav'];

%%
% Read mav file
% Data -> input signal
% fs -> sampling frequency

[Data,fs] = audioread(wav_file);

%% 
% Plot the fourier transform of the two channels from input signal
% Analysis window length in seconds
len = 0.040;
% Analysis window length in samples (power of 2 for faster FFT)
N0 = 2.^nextpow2(len*fs);
% Analysis window (even N and 'periodic' Hamming for constant overlap-add)
win = hamming(N0,'periodic');
stp = N0/2;
figure();
subplot(121);spectrogram(Data(:,1),hamming(N0),stp,N0,'yaxis');title([input_name, ' Left Channel of Input Signal'])
subplot(122);spectrogram(Data(:,2),hamming(N0),stp,N0,'yaxis');title([input_name, ' Right Channel of Input Signal'])

%%
% Process correlation and filter
y = Data;
for i = 1:N
    y = test(y, fs);
end

%%
audiowrite(background_output, y, fs);

%% 
%Plot the fourier transform of the two channels from output signal
figure();
subplot(121);spectrogram(y(:,1),hamming(N0),stp,N0,'yaxis');title([input_name, ' Left Channel of Output Signal'])
subplot(122);spectrogram(y(:,2),hamming(N0),stp,N0,'yaxis');title([input_name, ' Right Channel of Output Signal'])