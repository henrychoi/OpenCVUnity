function validate_CA(write_filter, run_exe, plot_fit)
if nargin < 1, write_filter = false; end
if nargin < 2, run_exe = true; end
if nargin < 3, plot_fit = true; end
%%
% My OSX CA input program saves the outputs to /tmp/ folder
M = 1024; % matched filter length
CAoutputDir = '/tmp/';
osxSrcDir = '~/github/OpenCVUnity/doc/CAInput/CH08_AUGraphInput/';
osxCAexe = strcat('~/DerivedData/', ...
    'CH08_AUGraphInput-fgzwyeogjkeirscegffwspzhghtm/Build/Products/',...
    'Release/', 'CH08_AUGraphInput');

unityAssetsDir = '~/github/OpenCVUnity/Assets/';
soundDir = strcat(unityAssetsDir, 'CubeShipsFree/Sound/');
streamingAssetDir = strcat(unityAssetsDir, 'StreamingAssets/');

unityDocDir = '~/github/OpenCVUnity/doc/';
fname = strcat(unityDocDir, 'FFT_pulse_conj.mat');
load(fname);

FFT_pulse_conj_vDSPz = FFT_pulse_conj(1:M); % Just take half
FFT_pulse_conj_vDSPz(1) = FFT_pulse_conj(1) + 1j*FFT_pulse_conj(M+1);
if write_filter
    %save(strcat(unityDocDir, 'FFT_pulse_conj.mat'), 'FFT_pulse_conj');
    % Using jsonlab-1.5, downloaded from Matlab file exchange
    savejson('', FFT_pulse_conj_vDSPz, ...
        strcat(streamingAssetDir, 'FFT_pulse_conj_vDSPz.json'));

    % Write the flattened form (row vector) of FFT_pulse_conj_vDSPz to the header
    FFT_pulse_conj_vDSPz_flattened = [ ...
        real(FFT_pulse_conj_vDSPz) imag(FFT_pulse_conj_vDSPz) ...
        ];
    csvwrite(strcat(osxSrcDir, 'FFT_pulse_conj_vDSPz_flattened.h'), ...
        FFT_pulse_conj_vDSPz_flattened);
end

% Now I should have the matched filter: FFT_pulse_conj (2M long)
temp = FFT_pulse_conj.';% .' = just transpose (rather than conj transpose)
if temp(end) ~= FFT_pulse_conj(end)
    error 'FFT_pulse_conj(end) wrong';
end
FFT_pulse_conj = temp; % do NOT flip the imag part!
clear temp;

% Play the sound and run the match filter program
[rchirp, Fs] = audioread('~/Music/rchirp.wav');
if run_exe
    %sound(rchirp, Fs);
    system(strjoin({osxCAexe, '--Nblock 12', ...
        '--Tcsv', strcat(CAoutputDir, 'CA_t.csv'), ...
        '--Xcsv', strcat(CAoutputDir, 'CA_x.csv'), ...
        '--Fcsv', strcat(CAoutputDir, 'CA_f.csv'), ...
        '--Ccsv', strcat(CAoutputDir, 'CA_c.csv'), ...
        '&'})); % background it!
    %java.lang.Thread.sleep(4)
    sound(rchirp, Fs);
end

if ~plot_fit
    return;
end
pause(1)
t = csvread(strcat(CAoutputDir, 'CA_t.csv'));
x = csvread(strcat(CAoutputDir, 'CA_x.csv'));
% block-by-block inv(FFT(x) * FFT(filter)) from CA
f = csvread(strcat(CAoutputDir, 'CA_f.csv'));% FFT(padded_x) .* FFT*(filter)
c = csvread(strcat(CAoutputDir, 'CA_c.csv'));% block-by-block correlation
N = length(x);

figure(1); clf;
subplot(311); plot(x); %semilogy(abs(x));

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
%FFT_CA = f(:,1) + 1j*f(:,2);

corr_CA = c(:,1);
d2_corr_CA = c(:,2);
med_corr_CA = c(:,3);
rel_corr_CA = c(:,4);

host_per_s  = (t(end,1) - t(1,1)) * (Fs / (2*M));
%world_per_s = (t(end,2) - t(1,2)) * (Fs / (2*M)); % all 0!

% What padded_x's half spectrum should be

% what block-by-block fast correlation result should be
corr_true = zeros(N,1);
overlap_save = zeros(M,1);
FFT_true = zeros(N,1);

for i=0:M:(N-M)
    rx_pad = [zeros(M,1); x((1:M) + i)];

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
    
    corr_true((1:M) + i) = overlap_save + fcorr(1:M);
    overlap_save = fcorr((1:M) + M); % save the partial for next round
end

ha1(1) = subplot(411); plot(x);

% vDSP forward FFT is 2x larger than correct values
 % vDSP inverse FFT is 2Mx larger than correct values (the filter length is
 % 2M)
 % In fast correlate (using forward and inverse FFT), the scale factor
 % would be 2*2*(2M) if vDSP is doing a forward FFT of the filter taps.
 % But I generate the FFT(filter) in Matlab, so the extra 2x for the filter
 % is is unnecessary.
corr_true = 2*(2*M) * real(corr_true);%corr should already be real
ha1(2) = subplot(412); % Show correlation
plot(corr_true((M+1):end));   hold on;
plot(corr_CA);
plot(corr_true((M+1):end) - corr_CA);
hold off;
legend('true correlation', 'CA correlation', 'correlation error')

med_corr_true1 = medfilt1(corr_true, 5);
corr_d1 = [corr_true(1); corr_true(1:end-1)];
corr_d2 = [corr_d1(1); corr_d1(1:end-1)];
corr_d3 = [corr_d2(1); corr_d2(1:end-1)];
corr_d4 = [corr_d3(1); corr_d3(1:end-1)];

A = max(corr_true, corr_d1);  B = min(corr_true, corr_d1);
C = max(corr_d2, corr_d3); D = min(corr_d2, corr_d3);
E = max(B, D);   D = min(B, D);
B = max(C, corr_d4);  C = min(C, corr_d4);
B = min(A, B);
A = max(E, C);   C = min(E, C);
E = max(A, B);   B = min(A, B);
D = max(D, C);
B = max(B, D);

rel_cor_d2_true = corr_d2 - B;

ha1(3) = subplot(413); % Show median of correlation
%plot(med_corr_CA);  hold on; plot(B((M+1):end));
%plot(med_corr_true1(((M+1):end) - 2) - med_corr_CA(:,1));
plot(rel_corr_CA); hold on;
plot(rel_cor_d2_true((M+1):end) - rel_corr_CA); hold on;
%plot(rel_cor_d2_true(((M+1):end) - 2)); hold on;
% plot(A((M+1):end));  hold on; plot(B((M+1):end));
% plot(med_corr_CA);
%plot([E((M+1):end) B((M+1):end)] - med_corr_CA);
hold off

d2Corr = [0; 0; diff(diff(corr_true))];

ha1(4) = subplot(414); % Show median of correlation
%plot(d2Corr((M+1):end)); hold on;
plot(d2_corr_CA); hold on;
plot(d2Corr((M+1):end) - d2_corr_CA);
hold off
legend('correlation double diff', 'double diff error');

% subplot(312); 
% %semilogy(abs(FFT_true));
% semilogy(abs(FFT_CA)); hold on;
% semilogy(abs(FFT_true - FFT_CA));
% hold off;
linkaxes(ha1, 'x');
