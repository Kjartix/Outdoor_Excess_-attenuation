using Plots, FFTW, DSP
plotly()

# Define chirp parameters
f_start = 16  # Starting frequency in Hz
f_end = 8000  # Ending frequency in Hz
T = 10         # Duration in seconds
Fs = 48000    # Sampling frequency in Hz
c = 343       # Speed of sound in m/s
d1 = 100      # Propagation distance 1 (m)
d2 = 100.3    # Propagation distance 2 (m)
octave_frequencies = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000]  # Hz

# Ensure consistent time vector
t = 0:1/Fs:T
t = t[1:end-1]  # Time vector
phase = f_start .* t .+ (f_end - f_start) .* (t.^2 /(2*T))
freqs = f_start .+ (f_end - f_start) .* (t / T)  # Linear frequency sweep
signal = sin.(2π .* phase .* t)

# Compute delayed signals
time_delay_1 = d1 / c  
time_delay_2 = d2 / c  

signal_d1 = sin.(2π .* phase .* t)
n = round(Int, time_delay_1 * Fs)  # 
signal_d2 = vcat(zeros(n), signal_d1)[1:length(signal_d1)]

x = signal_d1 .+ signal_d2  

# Compute spectrum
S = fft(signal)
X = fft(x)

# Compute spectrum
N = length(t)  # Number of points
positive_indices = 1:N÷2  # Take only positive frequencies
frequencies = fftfreq(N, Fs)[positive_indices]  # Extract positive frequencies

# Compute transfer function H = X / S
H_f = X ./ S
H_mag = abs.(H_f[positive_indices])
H_mag_dB = 20 * log10.(H_mag .+ eps())  # Convert to dB

# Create plots for time-domain signals
p1 = plot(t, x, label="Sum of signals at 100m", xlabel="Time (s)", ylabel="Amplitude", title="Sum of signals at 100m")
p2 = plot(t, signal_d1, label="Signal at 100m", xlabel="Time (s)", ylabel="Amplitude", title="Chirp Signal at 100m")
p3 = plot(t, signal_d2, label="Signal at 100.3m", xlabel="Time (s)", ylabel="Amplitude", title="Chirp Signal at 100.3m")

# Plot only positive frequencies
p4 = plot(frequencies, H_mag_dB, xlabel="Frequency (Hz)", ylabel="Magnitude (dB)", 
          title="Transfer Function Magnitude (dB)", xscale=:log10, label="Magnitude Response",
          xticks=(octave_frequencies, string.(octave_frequencies)))  # Custom x-axis ticks

# Add vertical lines at octave band frequencies
for f in octave_frequencies
    vline!([f], linestyle=:dash, color=:red, label="")
end

# Combine subplots
plot(p1, p2, p3, p4, layout=(4,1), size=(1100, 1100))
