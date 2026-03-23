% =========================================================================
% TRIDENT: 純流體物理縮放引擎 (Pure Hydro-Acoustic Scaling Engine)
% 零 k-Wave 依賴，專注於 Strouhal 頻移與 V^6 偶極子能量縮放
% =========================================================================
clear; clc; close all;

fprintf('啟動物理縮放引擎測試...\n');

%% 1. 定義經驗基準數據庫 (VA-111 Shkval 物理基準)
% 這些數據是我們縮放的「絕對錨點」
Baseline.V_ref = 100;       % 基準航速 (m/s)
Baseline.D_ref = 0.1;       % 基準鼻錐直徑 (m)
Baseline.SPL_ref = 195;     % 基準總輻射能量峰值 (dB)
Baseline.f_peak_ref = 500;  % 能量集中低頻峰值 (Hz)
Baseline.Power_Law = 60;    % 偶極子 V^6 定律 (60 * log10(V))

%% 2. 基礎訊號參數
fs = 500e3;           % 取樣率 500 kHz (足夠涵蓋 20kHz 尋標器頻帶)
duration = 0.1;       % 生成 0.1 秒的訊號來做頻譜分析
t_array = 0:(1/fs):(duration - 1/fs);
N_samples = length(t_array);
D_target = 0.1;       % 鼻錐直徑保持 10 cm

%% 3. 測試兩種不同航速的縮放結果
V_sprint = 100; % 衝刺航速 100 m/s
V_search = 60;  % 搜索航速 60 m/s

fprintf('生成 %d m/s 噪聲訊號...\n', V_sprint);
sig_100 = generate_scaled_noise(N_samples, fs, V_sprint, D_target, Baseline);

fprintf('生成 %d m/s 噪聲訊號...\n', V_search);
sig_60 = generate_scaled_noise(N_samples, fs, V_search, D_target, Baseline);

%% 4. 頻譜分析與對比繪圖 (證明縮放定律生效)
window = hanning(N_samples);
[Pxx_100, F] = pwelch(sig_100, window, [], N_samples, fs, 'onesided');
[Pxx_60, ~]  = pwelch(sig_60, window, [], N_samples, fs, 'onesided');

% 轉成 dB (基準值 1 uPa)
p_ref = 1e-6; 
Pxx_100_dB = 10 * log10(Pxx_100 / (p_ref^2));
Pxx_60_dB  = 10 * log10(Pxx_60 / (p_ref^2));

% 建立對比圖表
figure('Color', 'w', 'Position', [100, 100, 900, 500]);
semilogx(F, Pxx_100_dB, 'r', 'LineWidth', 2, 'DisplayName', '100 m/s');
hold on; grid on;
semilogx(F, Pxx_60_dB, 'b', 'LineWidth', 2, 'DisplayName', '60 m/s');

% 標示 20kHz 尋標器工作頻率的底噪差異
freq_check = 20000; 
[~, idx_20k] = min(abs(F - freq_check));
plot(freq_check, Pxx_100_dB(idx_20k), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility','off');
plot(freq_check, Pxx_60_dB(idx_20k), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility','off');

text(freq_check * 1.1, Pxx_100_dB(idx_20k) + 5, sprintf('%.1f dB', Pxx_100_dB(idx_20k)), 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');
text(freq_check * 1.1, Pxx_60_dB(idx_20k) - 5, sprintf('%.1f dB', Pxx_60_dB(idx_20k)), 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');

xlim([100, 100000]); 
ylim([100, 210]);
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power Spectral Density (dB re 1\muPa^2/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Hydro-Acoustic Scaling Engine Verification (Dipole V^6 Law)', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);

fprintf('縮放測試完成。請檢查繪出的頻譜對比圖。\n');

%% ========================================================================
% 核心縮放引擎函數 (Pure Scaling Logic)
% =========================================================================
function sig_noise = generate_scaled_noise(N_samples, fs, V_target, D_target, Baseline)
    % 1. 建立頻率軸
    df = fs / N_samples;
    f_full = (0:N_samples-1)' * df;
    
    f_interp = f_full;
    idx_neg = f_full > fs/2;
    f_interp(idx_neg) = fs - f_full(idx_neg); 
    
    % 2. 執行物理縮放 (Physics Scaling)
    % Strouhal 頻率偏移
    f_peak_current = Baseline.f_peak_ref * (V_target / Baseline.V_ref) * (Baseline.D_ref / D_target);
    % V^6 偶極子能量縮放
    SPL_peak_current = Baseline.SPL_ref + Baseline.Power_Law * log10(V_target / Baseline.V_ref);
    
    % 3. 建構目標頻譜包絡 (Envelope)
    target_curve_dB = zeros(size(f_interp));
    for i = 1:length(f_interp)
        f = f_interp(i);
        if f <= f_peak_current
            % 低頻段：能量快速爬升至峰值
            target_curve_dB(i) = SPL_peak_current - 20 * log10(f_peak_current / max(f, 1));
        else
            % 高頻段：嚴格套用 TBL 湍流衰減定律 (-16 dB / decade)
            target_curve_dB(i) = SPL_peak_current - 16 * log10(f / f_peak_current); 
        end
    end
    
    % 濾除直流與極高頻防呆
    target_curve_dB(1) = -200; 
    target_curve_dB(f_interp > fs/2.1) = -300; 
    
    % 4. 將 dB 轉回真實壓力幅值
    p_ref = 1e-6;
    P_linear_Pa2 = (p_ref^2) * 10.^(target_curve_dB / 10);
    Mag_Shape = sqrt(P_linear_Pa2); 
    
    % 5. 結合隨機相位生成時域訊號 (IFFT)
    rand_phase = 2 * pi * rand(N_samples, 1) - pi;
    half_len = floor(N_samples/2) + 1;
    Spec_Complex = zeros(N_samples, 1);
    Spec_Complex(1:half_len) = Mag_Shape(1:half_len) .* exp(1j * rand_phase(1:half_len));
    Spec_Complex(half_len+1:end) = conj(flipud(Spec_Complex(2:ceil(N_samples/2))));
    
    sig_raw = real(ifft(Spec_Complex));
    
    % 6. 精準能量標準化 (確保時域的變異數等於頻譜總能量)
    Total_Target_Power_Pa2 = sum(P_linear_Pa2(1:half_len)) * df; 
    Actual_Power = var(sig_raw);
    sig_noise = sig_raw * sqrt(Total_Target_Power_Pa2 / Actual_Power);
end