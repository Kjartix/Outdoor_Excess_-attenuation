using FFTW, Plots, LinearAlgebra, DSP

# Parameters
fs = 48000  # Sampling rate (Hz)
c = 343.0 
T=3  # Speed of sound (m/s)
d1 = 100.0  # Direct path distance (m)
d2 = 100.3  # Reflected path distance (m)
order = 4   # MLS order (length = 2^order - 1)
α = 1       # Reflection coefficient (attenuation factor)

# Compute time delays
time_delay_1 = d1 / c  
time_delay_2 = d2 / c  

# Generate MLS signal
function generate_mls(m::Int)
    N = 2^m - 1
    reg = ones(Int, m)  # Initial state
    seq = zeros(Int, N)

    for i in 1:N
        seq[i] = reg[end]  # Output bit
        feedback = xor(reg[1], reg[end])  # XOR 
        reg = [feedback; reg[1:end-1]]  # Shift register
    end

    return 2 .* seq .- 1  # Converting to -1 or 1
end

mls_signal = generate_mls(order)

# Compute frequency bins
freqs = fftfreq(length(mls_signal), 1/fs)

# Compute FFT of MLS signal
X = fft(mls_signal)

# Apply phase shifts for direct and reflected paths
phase_shift_1 = exp.(-im * 2 * π * freqs * time_delay_1)  # Direct path
phase_shift_2 = exp.(-im * 2 * π * freqs * time_delay_2)  # Reflected path

# Compute total received signal in frequency domain
Y = X .* phase_shift_1 + α * X .* phase_shift_2  # Sum of direct + reflection

# Compute inverse FFT to get the time-domain signal
mls_received = real(ifft(Y))

# Compute magnitude spectra
mls_spectrum = 20*log10.(abs.(X))  # Original MLS spectrum
received_spectrum = 20*log10.( abs.(Y))  # Received MLS spectrum

# Create subplots
plot(
    plot(mls_signal, label="Original MLS Signal", marker=:circle, title="Time Domain MLS Signal", xlabel="Sample Index", ylabel="Amplitude"),
    plot(mls_received, label="MLS with Reflection (100m + 100.3m)", linestyle=:dash, color=:red, title="Time Domain MLS with Reflection", xlabel="Sample Index", ylabel="Amplitude"),
    plot(freqs[1:div(end,2)], mls_spectrum[1:div(end,2)], label="Original MLS Spectrum", xlabel="Frequency (Hz)", ylabel="Magnitude", title="Frequency Spectrum of MLS Signal"),
    plot(freqs[1:div(end,2)], received_spectrum[1:div(end,2)], label="Received MLS Spectrum", linestyle=:dash, color=:red, xlabel="Frequency (Hz)", ylabel="Magnitude", title="Frequency Spectrum of MLS with Reflection"),
    layout=(2,2),
    size = (1200,1200)  # 2x2 grid of subplots
)

#sjekk spektrum på MLS
#Lenger enn impulsrespons, evt. ganske kort
#Sjekk litteratur