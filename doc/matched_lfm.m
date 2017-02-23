%%
function matched_lfm
COMPLEX_SIGNAL = false
ROUND_TRIP = false

if ROUND_TRIP, attenuation = 4; else attenuation = 2; end;

Fs = 48000;%44100;
Ts = 1/Fs;
f0 = 200; % lowest frequency of the chirp
b = 20000; % 20 kHz BW (use all available BW)
n = 2^13; % # samples to describe the original chirp
Nsample = 2^13; %Total samples I need to collect
c = 335; % speed of sound
taup = n * Ts;
Rmax = Nsample * Ts * c;
freqlimit = 0.5 * Fs;

c/b % range resolution

mu = b/taup;
scat_range = [2 5 10]; %[3.9 4 10];
scat_rcs = [1 1 1]; %[1 1.5 2];
nscat = length(scat_rcs);
winid = 0;
%eps = 1.0e-16;

s_ambient = audioread('singing.wav', (5 + [0 1]) * Nsample)';
RMS_ambient = sqrt(mean(s_ambient .^ 2));
A_u = 3 * RMS_ambient; %amplitude(s_u) = RMS(ambient)
sigma_r = 0.1 * RMS_ambient;

% time bandwidth product
time_B_product = b * taup;
if(time_B_product < 5 )
    fprintf('************ Time Bandwidth product is TOO SMALL ***************')
    fprintf('\n Change b and or taup')
  return
end
% check to ensure that scatterers are within recieve window
index = find(scat_range > Rmax);
if (index ~= 0)
    error('Error. Receive window is too large; or scatterers fall outside window')
end

% initialize input, output and replica vectors
x(nscat,1:Nsample) = 0.;
s_r = sigma_r * randn(1,Nsample) + s_ambient(1:Nsample);
replica(1:Nsample) = 0.;

% determine proper window
switch winid
    case 1
        win = hamming(Nsample)';
    case 2
        win = kaiser(Nsample, pi)';
    case 3
        win = chebwin(Nsample, 60)';
    otherwise
        win(1:Nsample) = 1.;
end

t = Ts * (0:Nsample-1);
t_xmit = t(1:n);

% Form the transmitted signal
if COMPLEX_SIGNAL
    replica(1:n) = A_u * exp(1j*2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
else %only the real
    replica(1:n) = A_u * cos(2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
end
freq = linspace(-freqlimit,freqlimit, n);
replica_mag = abs(fft(replica(1:n)));
FFT_replica = fft(replica, Nsample);

figure(1)
subplot(2,1,1); plot(t,real(replica)); grid
ylabel('real(replica)')
xlabel('t [s]')
subplot(2,1,2); plot(freq, fftshift(replica_mag)); grid
ylabel('|R|')
xlabel('[Hz]')

for i=1:50, sound(replica + s_r, Fs); end;

for j = 1:1:nscat
    range = scat_range(j);
    if ROUND_TRIP, range = 2*range; end;
    
    t_R = range/c; % nominal return time
    % x(:,k) is discrete sample of the returned signal, at time Ts*k.
    % Mahafza Equation (6.35): for t >= t_R
    % s_r ~ exp(2pi j fd*(t - t_R)) * exp(-2pi j f_0 t) * s_t(t - t_R)
    % where s_t is the original transmitted signal
    % Time period when the target will be subject to the signal I sent
    t_valid = (t >= t_R) & (t < (t_R + taup));
    t_xmit = t(t_valid) - t_R; % origin shifted by the return time
    if COMPLEX_SIGNAL
        x(j, t_valid) = A_u * exp(1j*2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
    else
        x(j, t_valid) = A_u * cos(2*pi * (f0 .* t_xmit + 0.5*mu .* t_xmit.^2));
    end
    
    % scale reflected signal and accumulate
    s_r = (scat_rcs(j) / (range)^attenuation) * x(j,:) + s_r;
end
s_r = s_r .* win; % Apply window

figure(2); semilogy(t, abs(replica) .^2, '--', t,abs(s_r).^2); %grid
title('Sound amplitude'); legend('xmit', 'received'); xlabel('[s]');

%
FFT_y = fft(s_r, Nsample);
correlation = ifft(FFT_y .* conj(FFT_replica)) / Nsample;

%correlation = xcorr(y, replica);
% For 2 vectors of length N, xcorr returns 2N-1 (all posible number of
% delays).  Negative delay means that the replica happened AFTER the
% output, which is acausal and therefore should be ignored.
%correlation = correlation(Nsample:end) / Nsample
figure(3)
distance = t*c;
if ROUND_TRIP, distance = distance/2; end;

semilogy(distance, abs(correlation))
%plot(t, abs(correlation))
xlabel ('Target relative distance [m]')
ylabel ('Xcorr')
grid

