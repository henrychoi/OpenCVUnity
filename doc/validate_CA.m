%%
% My OSX CA input program saves the outputs to /tmp/ folder
Fs = 44100;
M = 1024; % matched filter length
CAoutputDir = '/tmp/';
unityDocDir = '~/github/OpenCVUnity/doc/'
fname = strcat(unityDocDir, 'FFT_pulse_conj.mat');
load(fname);
% Now I should have the matched filter: FFT_pulse_conj (2M long)
temp = FFT_pulse_conj.';
if temp(end) ~= FFT_pulse_conj(end)
    error 'FFT_pulse_conj(end) wrong';
end
FFT_pulse_conj = temp; % do NOT flip the imag part!
clear temp;

t = csvread(strcat(CAoutputDir, 'CA_t.csv'));
x = csvread(strcat(CAoutputDir, 'CA_x.csv'));
% block-by-block inv(FFT(x) * FFT(filter)) from CA
f = csvread(strcat(CAoutputDir, 'CA_f.csv'));% Filter used by CA
c = csvread(strcat(CAoutputDir, 'CA_c.csv'));% block-by-block correlation
N = 0.5 * length(x);

figure(1); clf;
subplot(311); plot(x);

if false % Check matched filter
    %FFT_pulse_conj_CA = zeros(size(FFT_pulse_conj));
    FFT_pulse_conj_CA = f(:,1) + 1j*f(:,2);
    % This is how vDSP packs a real FFT; match it
    FFT_pulse_conj_vDSPz = FFT_pulse_conj(1:M); % Just take half
    FFT_pulse_conj_vDSPz(1) = FFT_pulse_conj(1) + 1j*FFT_pulse_conj(M+1);
    Filter_err = FFT_pulse_conj_vDSPz - FFT_pulse_conj_CA;
    subplot(212);
    %plot(imag(FFT_pulse_conj_vDSPz)); hold on;
    %plot(imag(FFT_pulse_conj_CA));
    plot(abs(Filter_err));
hold off;
end

% vDSP does not divide by 2, to save time, so do it here to compare against
% ground truth.
f = f/2;
FFT_CA = f(:,1) + 1j*f(:,2);

 % vDSP forward FFT is 2x larger than correct values
 % vDSP inverse FFT is 2Mx larger than correct values (the filter length is
 % 2M)
 % In fast correlate (using forward and inverse FFT), the scale factor
 % would be 2*2*(2M) if vDSP is doing a forward FFT of the filter taps.
 % But I generate the FFT(filter) in Matlab, so the extra 2x for the filter
 % is is unnecessary.
corr_CA = c / (2*(2*M));

host_per_s  = (t(end,1) - t(1,1)) * (Fs / (2*M));
%world_per_s = (t(end,2) - t(1,2)) * (Fs / (2*M)); % all 0!

% What padded_x's half spectrum should be

% what block-by-block fast correlation result should be
corr_true = zeros(2*N,1);
FFT_true = zeros(N,1);

for i=0:M:(N-M)
    rx_pad = x((1:2*M) + 2*i); % 0 padded x length: 2M

    % fast correlate: iFFT(fft(rx_pad) .* FFT_pulse_conj)
    fft_pad = fft(rx_pad);
    fft_mul = fft_pad .* FFT_pulse_conj;

    if true
        FFT_true((1:M) + i) = fft_mul(1:M);
        % vDSP moves the FFT @ Nyquist frequency to the imag(FFT @ 0)
        % to save space.  Match that behavior here.
        FFT_true(1+i) = FFT_true(1+i) + 1j * fft_mul(1+M);
    end
    fcorr = ifft(fft_mul);
    
    corr_true((1:(2*M)) + 2*i) = fcorr;
end
corr_true = real(corr_true);%corr should already be real

subplot(312); 
%semilogy(abs(FFT_true));
semilogy(abs(FFT_CA)); hold on;
semilogy(abs(FFT_true - FFT_CA));
hold off;

subplot(313);
%plot(corr_true);
plot(corr_CA);  hold on;
plot(corr_true - corr_CA);
hold off;

%corr_err = corr_true - corr_CA;
%subplot(313); plot(corr_err);

