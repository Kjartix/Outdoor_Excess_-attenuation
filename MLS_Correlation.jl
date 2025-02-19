using Random, FFTW, Plots, StatsBase, DSP
include("MLS.jl")
plotly()

# Generate white noise
d1 = 100.0  # Direct path distance (m)
d2 = 100.3 
fs = 96000  # Sampling frequency in Hz
T = 10.0  # Duration in seconds
N = Int(fs * T)  # Number of samples
n = round(Int, time_delay_1 * fs)  # 
time = range(0, stop=T, length=N)  # Time axis
order = 4
alpha = 1
c= 343 
time_delay_1 = d1 / c  
time_delay_2 = d2 / c  
lowcut = 16   # Lower cutoff frequency in Hz
highcut = 8000  
octave_frequencies = [16, 31.5, 63,125,250,500,1000,2000,4000,8000]  # Hz
#Generating white Noise
white_noise = randn(N)  # Generate white Gaussian noise

#Bandpass Filtering:
# Upper cutoff frequency in Hz
filter = digitalfilter(Bandpass(lowcut, highcut, fs=fs), Butterworth(100))
filtered_noise = filt(filter, white_noise) #Filtered


#Generate a repeating MLS sequence--------------------------
MLS = generate_mls(order)


MLS_len = length(MLS)
MLS_repeated = repeat(MLS, ceil(Int, N / MLS_len))[1:N]  # Trim to fit N samples
MLS_rep_filtered = filt(filter,MLS_repeated)#gir det meining?

encoded_MLS = MLS_rep_filtered .* filtered_noise #Filtered noise AND filtered MLS

#Magnitude spectrum------------------------------------------------------------------------------
frequencies = fftfreq(N, fs)[1:Int(N/2)]  
magnitude_spectrum_noise = 20*log10.( abs.(fft(filtered_noise))[1:Int(N/2)])

MLS_spectrum = 20*log10.( abs.(fft(encoded_MLS))[1:Int(N/2)])
#Autocorrelation:-------------------------------------------------------------------------------------------
AC_nonfiltered = autocor(MLS_repeated .- mean(MLS_repeated), 0:N÷2)# NON filtered signal
AC_nonfiltered = AC_nonfiltered ./AC_nonfiltered[1]  # Normalize


#Calcylating delays------------------------------------------------------------------------------
sample_delay_1 = round(Int, time_delay_1 * fs)  # Delay in number of samples
sample_delay_2 = round(Int, time_delay_2 * fs)  # Delay in number of samples

#delayed_signal_1 = encoded_MLS[time .- time_delay_1]
#delayed_signal_2 = zeros(N)

delayed_signal_1 = encoded_MLS
delayed_signal_2 = circshift(encoded_MLS,sample_delay_2)
#n = round(Int, time_delay_1 * fs)  # 
signal_d2 = vcat(zeros(n), delayed_signal_1)[1:length(delayed_signal_1)]
#signal_d2 = circshift(delayed_signal_1, sample_delay_2)

#--------------------------------------------------------------------------------------------------------
# x = MLS_rep_filtered
# x_d1 = MLS_rep_filtered()

MLS_sum = delayed_signal_1 .+ signal_d2


# p1=plot(time,delayed_signal_1)
# p2=plot(time,delayed_signal_2)

#FFT of the sum of the signals and the orignial signal
# #-------------------------------------------------------------------------------------------
positive_indices = 1:N÷2  # Take only positive frequencies
frequencies = fftfreq(N, fs)[positive_indices]  # Extract positive frequencies


X = fft(MLS_sum)
S = fft(encoded_MLS)
H_f = X ./ S
H_mag = abs.(H_f[positive_indices])
H_mag_dB = 20 * log10.(H_mag .+ eps())  # Convert to dB
#plot(p1,p2, layout=(2,1))


#dutille



# Plotting-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plot1 = plot(frequencies, magnitude_spectrum_noise, xlabel="Frequency (Hz)", ylabel="Magnitude",title="Magnitude Spectrum of White Noise", label="", lw=2, xlims=(10, 8500),xscale =:log10)
plot2 = plot(time, delayed_signal_1, xlabel="Time (s)", ylabel="Amplitude", title="White noise", label="", lw=1)
plot3 = plot(AC_nonfiltered, xlabel="Lag", ylabel="Autocorrelation", title="Autocorrelation Sequence of MLS encoded", label="", lw=2)
plot4 = plot(time,encoded_MLS, xlabel = "Time (s)", title = "ENcoded MLS signal", )
plot5 = plot(frequencies, MLS_spectrum,xlabel = "frequencies",title = "Encoded MLS spectrum",xscale =:log10)
plot6 = plot(time, delayed_signal_2, xlabel="Time (s)", ylabel="Amplitude", title="MLS Signal 100.3 m", label="", lw=1)
plot7 = plot(frequencies, H_mag_dB, xlabel="Frequency (Hz)", ylabel="Magnitude (dB)", title="Transfer Function Magnitude (dB)", xscale=:log10, label="Magnitude Response", xticks=(octave_frequencies, string.(octave_frequencies)))
for f in octave_frequencies
    vline!([f], linestyle=:dash, color=:red, label="")
end

# PLotting
plot(plot1,plot2,plot3,plot4,plot5,plot6,plot7,  layout=(7,1),size =(1200,1200))

