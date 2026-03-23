% =========================================================================
% TRIDENT: 剛性底線測試 (Solid Titanium Baseline & Tuning Fork Effect)
% Goal: 純鈦合金、無解耦器、正確流體聲源、4ms 長時間激盪
% =========================================================================
clear; clc; close all;

fprintf('Initializing Solid Titanium Baseline (Worst-Case Flanking)...\n');

%% 1. 建立高解析度運算網格
dx = 0.005; dy = 0.005; % 5 mm 
Nx = 1200; Ny = 400;   
kgrid = kWaveGrid(Nx, dx, Ny, dy);
X = kgrid.x; Y = kgrid.y;

%% 2. 定義純剛性介質 (100% Solid Titanium, NO Decoupler)
medium.sound_speed = 1500 * ones(Nx, Ny); % 水
medium.density = 1000 * ones(Nx, Ny); 
x_front = -2.0; 

% 中央陣列平台
central_platform_mask = (X >= x_front - 0.03) & (X < x_front) & (abs(Y) <= 0.025); 

% 車體外沿與主體 (全部連成一體)
cavitator_edge_mask = (X >= x_front) & (X < x_front + 0.05) & (abs(Y) <= 0.05); 
cone_start = x_front + 0.05; cone_L = 0.85; % 直接連接，無空隙
cone_radius = 0.05 + (0.2 - 0.05) * (X - cone_start) / cone_L;
cone_mask = (X >= cone_start) & (X < cone_start + cone_L) & (abs(Y) <= cone_radius);
hull_start = cone_start + cone_L;
hull_mask = (X >= hull_start) & (X <= 2.0) & (abs(Y) <= 0.2);

% 融合成一塊純鈦金屬
metal_mask = central_platform_mask | cavitator_edge_mask | cone_mask | hull_mask;
medium.sound_speed(metal_mask) = 5000; 
medium.density(metal_mask) = 4500; 

% 超空泡
cavity_L = 4.5; cavity_D = 0.8;
cavity_xc = x_front - 0.05 + cavity_L/2; 
cavity_mask = ((X - cavity_xc)/(cavity_L/2)).^2 + (Y/(cavity_D/2)).^2 <= 1;
gas_mask = cavity_mask & ~metal_mask;
medium.sound_speed(gas_mask) = 340; medium.density(gas_mask) = 1.2;  

medium.sound_speed = smooth(kgrid, medium.sound_speed, true);
medium.density = smooth(kgrid, medium.density, true);

%% 3. 生成真實寬頻物理噪聲 (TBL Pure Physics)
% 【延長模擬時間至 4ms，觀察完整的結構共振】
t_end = 4e-3; 
kgrid.makeTime(medium.sound_speed);
kgrid.t_array = 0:kgrid.dt:t_end;
fs = 1 / kgrid.dt; 
N_samples = length(kgrid.t_array);

Baseline.V_ref = 100; Baseline.D_ref = 0.1; Baseline.SPL_ref = 195; 
Baseline.f_peak_ref = 500; Baseline.Power_Law = 60;
[sig_source_Pa, ~] = generate_TBL_physics_noise(N_samples, fs, 100, 0.1, Baseline);

max_pressure_Pa = max(abs(sig_source_Pa));
sig_safe = sig_source_Pa / max_pressure_Pa; 

% 聲源精準放置於金屬外緣的水網格中
[~, edge_front_idx] = min(abs(kgrid.x_vec - x_front)); 
source_water_idx = edge_front_idx - 1; 
nose_y_idx = round(Ny/2);

source.p_mask = zeros(Nx, Ny);
source.p_mask(source_water_idx, nose_y_idx) = 1;       
source.p_mask(source_water_idx - 1, nose_y_idx) = 1;   
source.p = [sig_safe'; -sig_safe']; 

%% 4. 雙路徑模擬 (防止記憶體爆炸)
input_args = {'PlotSim', false, 'DataCast', 'single'};
[~, sensor_platform_idx] = min(abs(kgrid.x_vec - (x_front - 0.015))); % 平台中心

fprintf('\n--- PASS 1: Calculating Global SPL Field (4ms, ~6 mins) ---\n');
sensor1.mask = ones(Nx, Ny); 
sensor1.record = {'p_max'}; 
data_pass1 = kspaceFirstOrder2D(kgrid, medium, source, sensor1, input_args{:});

fprintf('\n--- PASS 2: Extracting Seeker Time-Series (4ms, ~6 mins) ---\n');
sensor2.mask = zeros(Nx, Ny);
sensor2.mask(sensor_platform_idx, nose_y_idx) = 1; 
sensor2.record = {'p'}; 
data_pass2 = kspaceFirstOrder2D(kgrid, medium, source, sensor2, input_args{:});

%% 5. 後處理與頻譜計算
p_ref = 1e-6; 

% 空間場
p_max_field_Pa = reshape(data_pass1.p_max, Nx, Ny) * max_pressure_Pa;
p_max_field_Pa(p_max_field_Pa == 0) = 1e-12;
SPL_field = 20 * log10(p_max_field_Pa / p_ref);

% 時間序列與 PSD
sig_received_Pa = data_pass2.p * max_pressure_Pa;

window = hanning(N_samples);
[Pxx_src, F_src] = pwelch(sig_source_Pa, window, [], N_samples, fs, 'onesided');
[Pxx_rec, F_rec] = pwelch(sig_received_Pa, window, [], N_samples, fs, 'onesided');

PSD_src_dB = 10 * log10(Pxx_src / (p_ref^2));
PSD_rec_dB = 10 * log10(Pxx_rec / (p_ref^2));

%% 6. 綜合繪圖
figure('Color', 'w', 'Position', [50, 50, 1400, 600]);

% --- 子圖 1: 空間 SPL 熱像圖 ---
subplot(1, 2, 1);
imagesc(kgrid.x_vec, kgrid.y_vec, SPL_field');
colormap(jet); caxis([110, max(SPL_field(:))]);
c = colorbar; ylabel(c, 'Max SPL (dB re 1\muPa)', 'FontWeight', 'bold');
hold on;
contour(kgrid.x_vec, kgrid.y_vec, cavity_mask', [0.5 0.5], 'w--', 'LineWidth', 1.5);
contour(kgrid.x_vec, kgrid.y_vec, metal_mask', [0.5 0.5], 'k-', 'LineWidth', 1.5);
plot(kgrid.x_vec(sensor_platform_idx), 0, 'wx', 'MarkerSize', 12, 'LineWidth', 2); 
xlim([-2.3, -1.8]); ylim([-0.2, 0.2]);
title('Solid Titanium Resonance (4ms Envelope)', 'FontSize', 14);
xlabel('Axial Distance (m)'); ylabel('Radial Distance (m)');
set(gca, 'YDir', 'normal'); 

% --- 子圖 2: 頻譜對比 ---
subplot(1, 2, 2);
semilogx(F_src, PSD_src_dB, 'r', 'LineWidth', 2, 'DisplayName', 'Source TBL Noise');
hold on; grid on;
semilogx(F_rec, PSD_rec_dB, 'b', 'LineWidth', 2, 'DisplayName', 'Received Noise (Solid Titanium Flanking)');

freq_sonar = 20000;
[~, idx_s] = min(abs(F_src - freq_sonar));
[~, idx_r] = min(abs(F_rec - freq_sonar));
plot(freq_sonar, PSD_src_dB(idx_s), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot(freq_sonar, PSD_rec_dB(idx_r), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
text(freq_sonar * 1.1, PSD_src_dB(idx_s), sprintf('%.1f dB', PSD_src_dB(idx_s)), 'Color', 'r', 'FontSize', 11);
text(freq_sonar * 1.1, PSD_rec_dB(idx_r), sprintf('%.1f dB', PSD_rec_dB(idx_r)), 'Color', 'b', 'FontSize', 11);

xlim([1000, 100000]); 
ylim([100, 220]);
xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Power Spectral Density (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Flanking Noise Transmission Spectrum (PSD)', 'FontSize', 14);
legend('Location', 'northeast');

fprintf('\n======================================\n');
fprintf('Solid Baseline Simulation Complete!\n');
fprintf('======================================\n');

%% ========================================================================
function [sig_noise, f_interp] = generate_TBL_physics_noise(N_samples, fs, V_target, D_target, Baseline)
    rand_phase = 2 * pi * rand(N_samples, 1) - pi;
    df = fs / N_samples;
    f_full = (0:N_samples-1)' * df;
    f_interp = f_full;
    idx_neg = f_full > fs/2;
    f_interp(idx_neg) = fs - f_full(idx_neg); 
    
    f_peak_current = Baseline.f_peak_ref * (V_target / Baseline.V_ref) * (Baseline.D_ref / D_target);
    SPL_peak_current = Baseline.SPL_ref + Baseline.Power_Law * log10(V_target / Baseline.V_ref);
    
    target_curve_dB = zeros(size(f_interp));
    for i = 1:length(f_interp)
        f = f_interp(i);
        if f <= f_peak_current
            target_curve_dB(i) = SPL_peak_current - 20 * log10(f_peak_current / max(f, 1));
        else
            target_curve_dB(i) = SPL_peak_current - 16 * log10(f / f_peak_current); 
        end
    end
    
    target_curve_dB(1) = -200; 
    
    p_ref = 1e-6;
    P_linear_Pa2 = (p_ref^2) * 10.^(target_curve_dB / 10);
    Mag_Shape = sqrt(P_linear_Pa2); 
    
    half_len = floor(N_samples/2) + 1;
    Spec_Complex = zeros(N_samples, 1);
    Spec_Complex(1:half_len) = Mag_Shape(1:half_len) .* exp(1j * rand_phase(1:half_len));
    Spec_Complex(half_len+1:end) = conj(flipud(Spec_Complex(2:ceil(N_samples/2))));
    
    sig_raw = real(ifft(Spec_Complex));
    Total_Target_Power_Pa2 = sum(P_linear_Pa2(1:half_len)) * df; 
    Actual_Power = var(sig_raw);
    sig_noise = sig_raw * sqrt(Total_Target_Power_Pa2 / Actual_Power);
    f_interp = f_full(1:half_len);
end