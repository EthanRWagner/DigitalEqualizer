%% Main Script
% USER MANUAL
% (1) Run Section -> "General Information and Specifications"
% (2) Run Section -> "Local Functions"
% (3) Comment/uncomment function calls to get the desired eq/rev settings
% (4) Run Section -> "Main Script"
%--------------------------------------------------------------------------
file = Zara;
[y,~] = audioread(file); % original audio file

% does nothing to the audio playback
% p_zara_arr = equalize_and_reverb_wavefile(file, EQ_settings_Nothing,...
% Dk_delays_msec_Nothing, alphak_gains_Nothing, ZarNothing);

% distorts and reverbs nonideally
p_zara_arr = equalize_and_reverb_wavefile(file, EQ_settings_Awful,...
Dk_delays_msec_Awful, alphak_gains_Awful, ZarAwful);

% very nice equalization and reverberation for peak experience
% p_zara_arr = equalize_and_reverb_wavefile(file, EQ_settings_Awesome,...
%     Dk_delays_msec_Awesome, alphak_gains_Awesome, ZarAwsome); 

fprintf("\nPlaying equalized...")
soundsc(p_zara_arr, fsample_CD);
pause(length(p_zara_arr)/fsample_CD);

% fprintf("\nPlaying original...")
% soundsc(y, fsample_CD)
% pause(length(y)/fsample_CD);

fprintf("\nDone...")
%% General Information and Specifications
fsample_CD = 44.1e3;
bands = [62.5, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3] ./ fsample_CD;

% input wave files
Zara19 = "ZaraExcerpt.wav";
Zara39S = "Zarathustra2.wav";
Zara39M = "Zarathustra1.wav";
Zara = "Zarathustra.wav";

% output wave files
ZarAwful = "Zarawful.wav";
ZarAwsome = "Zarawsome.wav";
ZarNothing = "Zarnothing.wav";
ZarEcho = "ZarEcho.wav";

% ZarAwful Specs
EQ_settings_Awful = [+12, -6, -12, -12, -6, +12, -6, +12, -6];
Dk_delays_msec_Awful = [250 400 520 660 750 1220];
alphak_gains_Awful = [0.7 0.6 0.5 0.33 0.2 0.8];

% ZarNothing Specs
EQ_settings_Nothing = [0,0,0,0,0,0,0,0,0];
Dk_delays_msec_Nothing = [];
alphak_gains_Nothing = [];

% ZarEcho Specs
Dk_delays_msec_Echo = [4000 10000];
alphak_gains_Echo = [1 1];

% ZarAwesome Specs
EQ_settings_Awesome = [4, 2, 7, 5, 4, 2, 10, 10, 10];
Dk_delays_msec_Awesome = [8850];
alphak_gains_Awesome = [0.75];

% FIR by Freq Sampling and windowing, long filter but no phase distortion
% amplitude range -12 to +12 dB

%% Local Functions
function eq_filter_hn = equalize_wavefile(eq_settings)
%   output_file = equalize_wavefile(input_file, eq_settings);
%   equalize_wavefile develops an FIR filter by frequency sampling and
%   applies the filter to the desired .wav file to output a processed
%   output .wav file.
%   The input parameters are:
%       eq_settings = 9 values representing the decibel levels for the
%       bands at 62.5 Hz, 125 Hz, 250 Hz, 500 Hz, 1 KHz, 2 KHz, 4 KHz, 8
%       KHz, and 16 KHz.
%
%   The output parameters are:
%       eq_filter_hn = FIR impulse response of EQ filter
%--------------------------------------------------------------------------

% Note: frequencies below 62.5 Hz and above 16 kHz take on magnitude at
% those frequencies
fsample_CD = 44.1e3;
bands = [62.5, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3] ./ fsample_CD;
% determine the frequency sample spacing --> provides M through 1/M
% relationship -> use CD fsample
deltaF = bands(1);
M = ceil(1/ deltaF);
if mod(M, 2) == 0
    M = M+1;
end
deltaF = (1/M);
eq_filter_hf = zeros(1, M);

% determine which samples are apart of the band center frequencies and add
% the respective decibel amplitude to those center frequencies
band_center_samps = zeros(1, length(bands));
for i=1:length(band_center_samps)
    band_center_samps(i) = round(bands(i)/deltaF);
    eq_filter_hf(band_center_samps(i)) = eq_settings(i);
end
% perform linear interpolation of dB values
% take the slope between center freq amplitudes and multiply sample number
% by coefficient and add to the left side center frequency amplitude
% linear interp = m*sample + lower_indexed_gain
% m = (db_high - db_low) / (samp_high - samp_low)
range = 1; % determines which set of center frequencies we are between
range_i = 0;
for i=1:(band_center_samps(end)-1)
    if (i >= band_center_samps(range))
        range=range+1;
        range_i = 0;
    end
    eq_filter_hf(i) = (eq_settings(range) - eq_settings(range-1)) ...
        / (band_center_samps(range) - band_center_samps(range-1)) *(range_i) + eq_settings(range-1);
    range_i = range_i + 1;
end
% make all frequencies higher than 16 KHz the same dB as 16 KHz
eq_filter_hf(band_center_samps(end):((M-1)/2) + 1) = eq_filter_hf(band_center_samps(end));
half_hn = eq_filter_hf(2:((M-1)/2 + 1));
eq_filter_hf = [eq_filter_hf(1) half_hn fliplr(half_hn)]; % make filter symmetric about Fs/2
% convert dB magnitudes into linear magnitudes
eq_filter_hf = 10.^(eq_filter_hf/20); % might use db2mag()

% convert frequency sampling array into an FIR filter unit sample response
[eq_filter_hn,~,~] = FIR_Filter_By_Freq_Sample(eq_filter_hf,1);

% apply Tukey window to reduce ringing and ripple
eq_filter_hn = real(eq_filter_hn.*tukeywin(length(eq_filter_hn)).');
[HF, Fd] = plot_DFT_mag(eq_filter_hn, fsample_CD, 1);
% fprintf("\nBefore Windowing");
% plot_freq_responses(Fd, eq_filter_hf, fsample_CD, 1);
% pause()
fprintf("\nAfter Windowing");
plot_freq_responses(Fd, HF, fsample_CD, 1);

end

function [echo_filter_hn]= echo_filter(Dk_delays_msec,alphak_gains,Fsample)
%   [echo_filter_hn]= echo_filter(Dk_delays_msec,alphak_gains,Fsample);
%   echo_filter constructs an FIR filter by having the original output
%   terms summated with any echo terms specified by the delay values and
%   gain coefficients in the input parameters.
%   The input parameters are:
%       Dk_delays_msec - the delay values (in ms) that an echo term should be
%           added
%       alphak_gains - the attenuation factor for the respective delay in
%           Dk_delays_msec
%       Fsample - sampling rate to determine the sample index from the delay
%           values
%
%   The output parameters are:
%       echo_filter_hn = the FIR filter impulse response to produce the
%           desired echo effect
%--------------------------------------------------------------------------

if(~isempty(Dk_delays_msec))
    % create array of zeros based on fsample and last Dk value
    echo_filter_hn = zeros(1, floorDiv(Dk_delays_msec(end)/1000*Fsample, 1));
    echo_filter_hn(1) = 1; % include x[n]
    % for each delay in Dk
    for i=1:length(Dk_delays_msec)
        % convert ms to s then to the sample index
        samp_num = floorDiv(Dk_delays_msec(i)/1000*Fsample, 1);
        % add fractional gain to hn at sample index above
        echo_filter_hn(samp_num) = alphak_gains(i);
    end
else
    echo_filter_hn = [1]; % return y[n]=x[n] if no delays
end
end

function Processed_wav = equalize_and_reverb_wavefile(inwavfilename, ...
    EQdBsettings, Dk_delays_msec, alphak_gains, outwavfilename)
%   Processed_wav = equalize_and_reverb_wavefile(inwavfilename,...
%       EQdBsettings, Dk_delays_msec, alphak_gains, outwavfilename);
%   equalize_and_reverb uses equalize_wavefile and echo_filter to equalize
%   and reverberate a wavefile to enhance (or diminish) its acoustic
%   properties
%   The input paramters are:
%        inwavfilename - the input file name
%        EQdBsettings - the db amplitudes at the 9 frequency bands described in
%           equalize_waveform
%        Dk_delays_msec - the delay values (in ms) that an echo term should be
%           added
%        alphak_gains - the attenuation factor for the respective delay in
%           Dk_delays_msec
%        outwavfilename - the output file name to be created
%
%   The output parameters are:
%       Processed_wav - the output .wav file
%--------------------------------------------------------------------------
% read in audio file
[y, fs] = audioread(inwavfilename);

% obtain equalizer and reverberation filters
eq_filter_hn = equalize_wavefile(EQdBsettings);
echo_filter_hn = echo_filter(Dk_delays_msec, alphak_gains, fs);

% use fftconv() to convolve both filters together
total_filter_hn = fftconv(eq_filter_hn, echo_filter_hn);

% determine if the input .wav file is single or 2-channel and process each
% channel individually if necessary
channels = 2; % default to stereo
dim_y = size(y);
if(dim_y(2) == 1)
    channels = 1;
end
% use fftconv() to apply resulting filter to the .wav file
for i=1:channels
    B = fftconv(y(:, i).', total_filter_hn).';
    y(1:length(B), i) = B;
end
% normalize .wav file with max(abs(input.wav))
max_val = max(abs(y(:)));
y = y ./ max_val;

% write .wav file and output processed data
audiowrite(outwavfilename, y, fs);
Processed_wav = y;

end

function yn = fftconv(xn, hn)
%function yn = fftconv( xn, hn )
%   where the input arguments are:
%       xn = time domain samples of a discrete-time input signal
%       hn = time domain samples of a D-T system Unit Sample Response 
%   and the output is:
%       yn = time domain samples of a discrete-time system output signal

LIMIT = 1000; % limit of samples for figure rendering

% compute zero padding length
x_length = length(xn); % length of x
n_x = 0:(x_length-1); % sample array for x
h_length = length(hn); % length of h
n_h = 0:(h_length-1); % sample array for h

% what the ACTUAL y length will be
y_length = x_length + h_length - 1;
n_y = 0:(y_length-1);

% length to zero pad sequences to a power of 2 to speed up fft
npow2_length = 2^nextpow2(y_length);

x_zs = zeros(npow2_length - x_length, 1).'; % num of 0s to append to x
h_zs = zeros(npow2_length - h_length, 1).'; % num of 0s to append to h

x_p = [xn x_zs]; % x zero padded sequence
h_p = [hn h_zs]; % h zero padded sequence

% take fft of each sequence
xk = fft(x_p); % dft of zero padded xn
n_x_k = (0:(length(xk)-1)) ./ length(xk); % sample values of xk
hk = fft(h_p); % dft of zero padded hn
n_h_k = (0:(length(hk)-1)) ./ length(hk); % sample values of hk

% complex multiply X(k).*H(k)
yk = xk.*hk; % freq domain multiplication
n_y_k = (0:(length(yk)-1)) ./ length(yk); % sample values of yk

% Y(k) is product and y[n] is ifft
yn_zp = ifft(yk); % inverse fft to get yn
yn = yn_zp(1:y_length);

end

function [hn,HF,F] = FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum)
%[hn,HF,F] = FIR_Filter_By_Freq_Sample(HF_mag_samples,figurenum);
%   The function creates a casual Type I filter by way of frequency sampling.
%   The input parameters are
%       HF_mag_samples – H[k] Magnitude response samples for desired filter
%       figurenum - Figure # to plot frequency responses
%   The output parameters are
%       hn - impulse response of filter (same length as HF_mag_samples)
%       HF - complex frequency response of filter
%           (estimated H(F) values found by FFT or freqz)
%       F – digital frequency values corresponding to the estimated H(F) 
%           values
%--------------------------------------------------------------------------
M = length(HF_mag_samples); % get filter length
zp_length = 2048; % zero padded length
k = 0:(M-1); % k values for iteration
Fx = linspace(0, (M-1)/M, M); % digital frequency value
phase = -pi*(M-1)/M*k; % obtain phase values via Type 1 phase equation
polar_phase = exp(i*phase); % e raised to the phase value

HFx = HF_mag_samples.*polar_phase; % mag * e^phase to get freq resp values
hn = ifft(HFx, M); % get time domain sequence of impulse response

zp = zeros(1, (zp_length - M)); % zero padded sequence to be appended
hn_zp = [hn zp]; % zero pad the impulse response
HF = fft(hn_zp); % get freq repsonse values for continuous plot
F = (0:(zp_length-1)) ./ zp_length; % get freq values

end

function plot_freq_responses(Fd, HF, fsample, figure_num)
% plot_freq_responses(Fd, HF, fsample, figure_num);
%   Detailed explanation goes here

figure(figure_num)
% Plot the Magnitude Response 
subplot(2,1,1)  % Display plots in 2 rows / 1 column; This is the 1st plot.

% Plot the magnitude of HF on a linear scale
plot(Fd, abs(HF))
grid on
xlabel('Digital Frequency  F (cycles/sample)')
ylabel('Magnitude Response')
title('Frequency Response of Filter')

% Plot the Phase Response below the Magnitude Response
subplot(2,1,2) % Display plots in 2 rows / 1 column; This is the 2nd plot.

% Plot the Phase Angle vs Frequency     
plot(Fd, angle(HF)/pi)        % Normalize angle radian values by pi radians
grid on
xlabel('Digital Frequency  F (cycles/sample)')
ylabel('Phase Response /pi')

figure(figure_num + 1)
% Plot the Magnitude Response 
subplot(2,1,1)
% Plot the magnitude of HF on a linear scale
plot(Fd, 20*log10(abs(HF)))
grid on
xlabel('Digital Frequency  F (cycles/sample)')
ylabel('Magnitude Response (dB)')
title('Frequency Response of Filter')

% Plot the Phase Response below the Magnitude Response
subplot(2,1,2) % Display plots in 2 rows / 1 column; This is the 2nd plot.

% Plot the Phase Angle vs Frequency     
plot(Fd, angle(HF)/pi)        % Normalize angle radian values by pi radians
grid on
xlabel('Digital Frequency  F (cycles/sample)')
ylabel('Phase Response /pi')

end

function [DFTx,Fd] = plot_DFT_mag(x,fsample, figure_num)
% function [DFTx,Fd] = plot_DFT_mag(x,fsample, figure_num);
%   where the input arguments are:
%       x = time domain samples of a discrete-time sequence
%       fsample = sampling frequency (samples / second)
%       figure_num = number of the figure to use for the two plots
%   and the outputs are:
%       DFTx = DFT spectrum values (complex)[same # samples as x]
%       Fd = Digital frequency values for each DFT sample

DFTx = fft(x); % DFT of x
num_of_samples = length(DFTx); % number of samples
Fd = (0:(num_of_samples - 1)) ./ num_of_samples; % digital freq value

% % plot magnitude responses
% figure(figure_num);
% subplot(2, 1, 1); % digital freq plot
% stem(Fd, abs(DFTx) / num_of_samples, "."); % stem plot with dot markers
% xlabel("Digital Frequency (cycles/sample)");
% ylabel("DFT Magnitude");
% title("DFT Magnitude Spectrum");
% 
% subplot(2, 1, 2); % analog freq plot
% stem(Fd.*fsample, abs(DFTx) / num_of_samples, "."); % stem plot with dot markers
% xlabel("Analog Frequency (Hz)");
% ylabel("DFT Magnitude");

end