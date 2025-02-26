%% Parameters
order = 7; % Order of MLS sequence
fs = 48000; % Sampling frequency (Hz)
duration = 10; % Duration of the signal in seconds
c = 343; % Speed of sound (m/s)

%% Bandpass filter parameters
low_cutoff = 16;  % Hz
high_cutoff = 8000; % Hz
filter_order = 5; % Filter order for Butterworth filter

%% Design Bandpass Filter
[b, a] = butter(filter_order, [low_cutoff high_cutoff] / (fs/2), 'bandpass');

%% Generate MLS sequence
mls_sequence = mls(order); % Generates an MLS sequence of length 2^order - 1

%% Calculate the required length to match the white noise sequence
required_length = fs * duration;

%% Repeat the MLS sequence to match the required length
repeated_mls = repmat(mls_sequence, ceil(required_length / length(mls_sequence)), 1);
repeated_mls = repeated_mls(1:required_length); % Truncate to match exact length

%% Generate white noise sequence
white_noise = randn(required_length, 1); % Gaussian white noise

%% Normalize both sequences to have the same RMS level
repeated_mls = repeated_mls / rms(repeated_mls);
white_noise = white_noise / rms(white_noise);

%% Modulate the MLS sequence into the white noise
modulated_signal = white_noise + repeated_mls;

%% Simulate propagation delays
d1 = 100; % Distance 1 in meters
d2 = 100.3; % Distance 2 in meters

t1 = d1 / c; % Delay in seconds for 100m
t2 = d2 / c; % Delay in seconds for 100.3m

n_delay1 = round(t1 * fs); % Convert delay to samples
n_delay2 = round(t2 * fs); % Convert delay to samples

%% Apply attenuation
A1 = 1 / d1;
A2 = 1 / d2;

%% Create delayed versions of the modulated signal
signal1 = A1 * [zeros(n_delay1,1); modulated_signal(1:end-n_delay1)];
signal2 = A2 * [zeros(n_delay2,1); modulated_signal(1:end-n_delay2)];

%% Sum the two received signals (Before filtering)
received_signal = signal1 + signal2;

%% Apply Bandpass Filter to the received signal
filtered_signal1 = filtfilt(b, a, signal1);
filtered_signal2 = filtfilt(b, a, signal2);
filtered_received_signal = filtered_signal1 + filtered_signal2;

%% Compute Welch’s Power Spectral Density Estimate
window_length = 4096;  % Choose window length
overlap = 2048;        % 50% overlap
nfft = 4096;           % FFT size

[Pxx, freq_welch] = pwelch(filtered_received_signal, hamming(window_length), overlap, nfft, fs);

%% Compute Autocorrelation of the MLS Sequence
mls_autocorr = xcorr(mls_sequence, 'unbiased');
mls_lags = (-length(mls_sequence)+1:length(mls_sequence)-1) / fs;

%% Time axis for plotting
t = (0:length(repeated_mls)-1) / fs;

%% 🎨 Plotting Results
figure;

% 1️⃣ Plot MLS Sequence
subplot(4,1,1);
stem(mls_sequence, 'k');
title('MLS Sequence');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% 2️⃣ Plot Outgoing Modulated Signal (MLS + White Noise)
subplot(4,1,2);
plot(t, modulated_signal, 'b');
title('Outgoing Signal (Modulated MLS)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% 3️⃣ Plot Autocorrelation of MLS
subplot(4,1,3);
plot(mls_lags, mls_autocorr, 'r');
title('Autocorrelation of MLS Sequence');
xlabel('Lag (s)');
ylabel('Autocorrelation');
grid on;

% 4️⃣ Plot Welch’s Power Spectral Density (PSD)
subplot(4,1,4);
semilogx(freq_welch, 10*log10(Pxx), 'm');
title('Power Spectral Density (Welch’s Method)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([10 fs/2]);
