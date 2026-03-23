% =========================================================================
% TRIDENT: Tactical Acoustic Spectrum (Pure Simulated Engine Output)
% Ultimate Calibrated Version (PSD Bin Width & Single-Sided Energy Fixed)
% =========================================================================
clear; clc; close all;

%% 1. 內部驅動參數 (Baseline & Tactical Scaling)
V1 = 10; D1 = 0.0254;        
skid_freq = [50, 100, 200, 500, 1000, 3000, 6000, 10000, 20000, 40000];
skid_ampl = [150, 170, 208, 185, 195, 212, 205, 180, 165, 150];

V2 = 100; D2 = 0.1;           
va111_freq = skid_freq .* (V2 / V1) .* (D1 / D2);
va111_ampl = skid_ampl + 40 * log10(V2 / V1);

%% 2. 呼叫 TRIDENT 隨機物理引擎 
fs = 200000;
T_sim = 1.0; 
N_samples = round(T_sim * fs);

fprintf('Generating mathematically calibrated stochastic noise field...\n');
[sig_tail, env_curve, f_env] = generate_spiky_noise(N_samples, fs, va111_freq, va111_ampl);

%% 3. 計算真實的功率譜密度 (PSD) - 修正 'psd' 參數
Nfft = 16384; 
window = hanning(Nfft);
% 【關鍵修正 1】：使用 'psd' 確保繪製的是能量密度，避免頻寬陷阱帶來的 11 dB 誤差
[Pxx_tail, F_psd] = pwelch(sig_tail, window, Nfft/2, Nfft, fs, 'psd'); 
SPL_tail_sim = 10*log10(Pxx_tail);

%% 4. 繪製高學術標準圖表 (只保留曲線與漸近線)
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
hold on; grid on;

% 繪製理論包絡線/漸近線 (黑色虛線)
plot(f_env, env_curve, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Theoretical $f^{-2}$ Envelope');

% 繪製引擎生成的真實噪聲頻譜 (藍色尖刺線)
plot(F_psd, SPL_tail_sim, 'b', 'LineWidth', 0.8, 'DisplayName', 'Simulated Tail Source PSD (Time-Domain Synthesis)');

% 圖表格式設定
set(gca, 'XScale', 'log');
xlim([10, 100000]);
ylim([120, max(va111_ampl) + 15]);
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
% 【嚴謹細節】：PSD 的標準單位是 dB re 1\muPa^2/Hz
ylabel('Sound Pressure Level (dB re 1\muPa^2/Hz)', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Simulated Supercavitating Tactical Noise Spectrum (100 m/s Regime)', 'FontSize', 14);

leg = legend('Location', 'northeast', 'FontSize', 11);
set(leg, 'Interpreter', 'latex');

% =========================================================================
% 隨機物理引擎 (包含連續譜與洛倫茲尖刺，單邊能量校準)
% =========================================================================
function [sig_noise, target_curve_dB_full, f_interp] = generate_spiky_noise(N_samples, fs, peak_freqs, peak_ampls)
    rand_phase = 2 * pi * rand(N_samples, 1) - pi;
    df = fs / N_samples;
    f_full = (0:N_samples-1)' * df;
    
    f_interp = f_full;
    idx_neg = f_full > fs/2;
    f_interp(idx_neg) = fs - f_full(idx_neg);
    
    % 基線連續譜
    base_level = max(peak_ampls) - 45; 
    target_curve_dB_full = base_level * ones(size(f_interp));
    
    % 高頻段慣性次區衰減 (-20dB/dec)
    f_decay_start = peak_freqs(end);
    idx_high = f_interp > f_decay_start;
    target_curve_dB_full(idx_high) = base_level - 20*log10(f_interp(idx_high) / f_decay_start);
    
    % 注入洛倫茲諧波尖刺
    for i = 1:length(peak_freqs)
        Q = 40; 
        bw = peak_freqs(i) / Q;
        spike_profile = peak_ampls(i) - 10*log10(1 + 4*((f_interp - peak_freqs(i))/bw).^2);
        target_curve_dB_full = max(target_curve_dB_full, spike_profile);
    end
    
    target_curve_dB_full(1) = -200; 
    
    % 轉換回時域
    P_linear = 10.^(target_curve_dB_full / 10);
    Mag_Shape = sqrt(P_linear); 
    Spec_Complex = Mag_Shape .* exp(1j * rand_phase);
    sig_raw = real(ifft(Spec_Complex));
    
    % 【關鍵修正 2】：能量標準化嚴格限制在單邊頻譜，消除 3dB 折疊誤差
    half_len = floor(N_samples/2) + 1;
    Total_Target_Power = sum(P_linear(1:half_len)) * df; 
    Actual_Power = var(sig_raw);
    
    if Actual_Power > 0
        sig_noise = sig_raw * sqrt(Total_Target_Power / Actual_Power);
    else
        sig_noise = zeros(size(sig_raw));
    end
    sig_noise = sig_noise(:);
    
    % 確保回傳的理論曲線對齊正確的頻率長度
    target_curve_dB_full = target_curve_dB_full(1:half_len);
    f_interp = f_full(1:half_len);
end