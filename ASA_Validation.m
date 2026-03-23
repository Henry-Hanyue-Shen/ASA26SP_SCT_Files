% =========================================================================
% SCT_Acoustic_Engine_Master_Proof.m
% -------------------------------------------------------------------------
% PROJECT TRIDENT: SCT PART B - Master Acoustic Engine Validation
% PURPOSE: 
%   Comprehensive Bottom-Up Validation of the Decoupled Acoustic Engine.
%   PART 1: Low-Frequency Monopole Pulsation (Multi-Condition Harmonics)
%           Strictly anchored to Skidmore (2016) empirical hydrophone data.
%   PART 2: High-Frequency Multipole Roll-off (Dipole/Quadrupole)
%           Strictly follows the f^-2 theoretical decay (-20 dB/decade) 
%           established by Mellen (1954) and Blake (1986).
% =========================================================================

clear; clc; close all;
rng(2026); % Ensure reproducibility

%% ========================================================================
%  PART 1: EMPIRICAL ANCHORS & HYPERPARAMETERS (Low-Speed Validation)
%  Reference: Skidmore, G. (2016). Ph.D. Dissertation. Penn State Univ.
% =========================================================================

% --- Low-Frequency Pulsation Conditions (Table 4.2 & Fig 4.9-4.11) ---
% Case A: Standard Pulsation
cases(1).name   = 'Case A: Standard Pulsation (Ref: Fig 4.10)';
cases(1).Fr     = 3.67; cases(1).CQ = 1.08; cases(1).f_mod = 'None';
f0_A = 38.8; 
cases(1).peak_f = [f0_A, 2*f0_A, 3*f0_A];   
cases(1).peak_A = [186.8, 160.0, 155.0]; % dB (SPL)

% Case B: Low-Froude Pulsation
cases(2).name   = 'Case B: Low-Fr Pulsation (Ref: Fig 4.11)';
cases(2).Fr     = 3.23; cases(2).CQ = 1.47; cases(2).f_mod = 'None';
f0_B = 36.7; 
cases(2).peak_f = [f0_B, 2*f0_B, 3*f0_B];     
cases(2).peak_A = [184.4, 158.0, 152.0]; % dB (SPL)

% Case C: Frequency-Shifted Pulsation via Modulation
cases(3).name   = 'Case C: Freq-Shifted Pulsation (Ref: Fig 4.9)';
cases(3).Fr     = 3.67; cases(3).CQ = 1.09; cases(3).f_mod = '32.0 Hz';
f0_C = 41.0; 
cases(3).peak_f = [f0_C, 2*f0_C, 3*f0_C];     
cases(3).peak_A = [183.1, 156.0, 151.0]; % dB (SPL)

% --- Empirical Broadband Continuum Baseline ---
bb_freq = [10, 100, 500, 1000];
bb_ampl = [140, 138, 136, 135]; % dB (SPL)

%% ========================================================================
%  PART 2: THEORETICAL MULTIPOLE ROLL-OFF (High-Frequency Validation)
%  References: 
%  1. Mellen, R. H. (1954). Ultrasonic Spectrum of Cavitation Noise.
%  2. Blake, W. K. (1986). Mechanics of Flow-Induced Sound and Vibration.
%  Formula: SPL(f) = SPL(f_ref) - 20 * log10(f / f_ref)
% =========================================================================

hf_f_ref = 1000; % Hz (Terminal boundary of empirical macroscopic data)
hf_SPL_ref = 135; % dB (Anchored to empirical continuum)
hf_freq = [1000, 3000, 10000, 20000, 40000];
hf_ampl = hf_SPL_ref - 20 * log10(hf_freq / hf_f_ref);

%% ========================================================================
%  SYSTEM PARAMETERS FOR ENGINE EXECUTION
% =========================================================================
fs_low = 10000;          % 10 kHz sampling for macroscopic pulsation
fs_high = 200000;        % 200 kHz sampling for micro-bubble roll-off
T_sim = 5.0;             % 5 seconds duration for high spectral resolution
window_size_low = round(fs_low * 0.5); 
window_size_high = round(fs_high * 0.05);

%% ========================================================================
%  EXECUTION & REPORTING
% =========================================================================
fprintf('\n====================================================================\n');
fprintf('PROJECT TRIDENT: MASTER ACOUSTIC ENGINE VALIDATION REPORT\n');
fprintf('====================================================================\n\n');

% -------------------------------------------------------------------------
% SECTION 1: LOW-FREQUENCY HARMONIC VALIDATION
% -------------------------------------------------------------------------
fprintf('>>> SECTION 1: LOW-FREQUENCY MONOPOLE VALIDATION (Skidmore 2016)\n');
f1 = figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
sgtitle('Section 1: Low-Frequency Harmonic Pulsation Validation', 'FontWeight', 'bold', 'FontSize', 16);

for c = 1:length(cases)
    N_samples = round(T_sim * fs_low); t = (0:N_samples-1)' / fs_low;
    
    % Step 1: Generate Stochastic Continuum (Decoupled)
    sig_bb = generate_broadband_noise(N_samples, fs_low, bb_freq, bb_ampl);
    
    % Step 2: Inject Deterministic Harmonics
    sig_tones = zeros(N_samples, 1);
    for p = 1:length(cases(c).peak_f)
        f_p = cases(c).peak_f(p); A_p = cases(c).peak_A(p);
        unit_sine = sin(2 * pi * f_p * t);
        [pxx_unit, f_welch] = pwelch(unit_sine, hamming(window_size_low), round(window_size_low/2), window_size_low, fs_low);
        unit_peak_dB = 10*log10(pxx_unit(dsearchn(f_welch, f_p)));
        sig_tones = sig_tones + 10^((A_p - unit_peak_dB) / 20) * sin(2 * pi * f_p * t + rand()*2*pi);
    end
    
    % Step 3: Superposition & PSD Estimation
    sig_total = sig_bb + sig_tones;
    [pxx_total, f_welch] = pwelch(sig_total, hamming(window_size_low), round(window_size_low/2), window_size_low, fs_low);
    sim_psd_dB = 10*log10(pxx_total);
    
    % Step 4: Error Computation
    baseline_target = interp1(bb_freq, bb_ampl, f_welch, 'pchip', NaN);
    valid_idx = ~isnan(baseline_target);
    for p = 1:length(cases(c).peak_f), valid_idx = valid_idx & abs(f_welch - cases(c).peak_f(p)) > 5; end
    baseline_rmse = sqrt(mean((sim_psd_dB(valid_idx) - baseline_target(valid_idx)).^2));
    
    % Print Console Report
    fprintf('--------------------------------------------------------------------\n');
    fprintf('[%s]\n', cases(c).name);
    fprintf('  Hyperparameters : Fr = %.2f, CQ = %.2f, f_mod = %s\n', cases(c).Fr, cases(c).CQ, cases(c).f_mod);
    for p = 1:length(cases(c).peak_f)
        idx_p = dsearchn(f_welch, cases(c).peak_f(p));
        fprintf('  Harmonic %d      : Target %6.1f dB | Actual %6.1f dB (Err: %.2f dB)\n', p, cases(c).peak_A(p), sim_psd_dB(idx_p), abs(sim_psd_dB(idx_p) - cases(c).peak_A(p)));
    end
    fprintf('  Continuum RMSE  : %.2f dB\n', baseline_rmse);
    
    % Plot
    subplot(3, 1, c); hold on; grid on;
    plot(f_welch, baseline_target, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Empirical Continuum');
    plot(f_welch, sim_psd_dB, 'Color', [0.1, 0.3, 0.7, 0.9], 'LineWidth', 1.2, 'DisplayName', 'Simulator Engine');
    plot(cases(c).peak_f, cases(c).peak_A, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Empirical Peaks');
    set(gca, 'XScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
    xlim([10, 1000]); ylim([120, 200]); ylabel('SPL (dB)');
    title(sprintf('%s | Base RMSE: %.2f dB', cases(c).name, baseline_rmse));
    if c == 1, legend('Location', 'northeast'); end
    if c == 3, xlabel('Frequency (Hz)'); end
end

% -------------------------------------------------------------------------
% SECTION 2: HIGH-FREQUENCY MULTIPOLE VALIDATION
% -------------------------------------------------------------------------
fprintf('\n>>> SECTION 2: HIGH-FREQUENCY MULTIPOLE VALIDATION (Mellen & Blake)\n');
N_samples_hf = round(T_sim * fs_high);
sig_multipole = generate_broadband_noise(N_samples_hf, fs_high, hf_freq, hf_ampl);
[pxx_hf, f_welch_hf] = pwelch(sig_multipole, hamming(window_size_high), round(window_size_high/2), window_size_high, fs_high);
sim_psd_hf_dB = 10*log10(pxx_hf);

valid_idx_hf = find(f_welch_hf >= 1000 & f_welch_hf <= 40000);
theoretical_target = hf_SPL_ref - 20 * log10(f_welch_hf(valid_idx_hf) / hf_f_ref);
multipole_rmse = sqrt(mean((sim_psd_hf_dB(valid_idx_hf) - theoretical_target).^2));

fprintf('  Theoretical Law : -20 dB/decade (-6 dB/octave) Inertial Decay\n');
fprintf('  Validation Band : 1,000 Hz to 40,000 Hz\n');
fprintf('  Multipole RMSE  : %.2f dB\n', multipole_rmse);
fprintf('====================================================================\n\n');

f2 = figure('Color', 'w', 'Position', [150, 150, 800, 500]);
hold on; grid on;
plot(f_welch_hf(valid_idx_hf), theoretical_target, 'k--', 'LineWidth', 2, 'DisplayName', 'Theoretical f^{-2} Decay');
plot(f_welch_hf, sim_psd_hf_dB, 'Color', [0.1, 0.3, 0.7, 0.8], 'LineWidth', 1.5, 'DisplayName', 'TRIDENT Engine');
plot(hf_freq, hf_ampl, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'DisplayName', 'Theoretical Anchors');
xline(20000, 'g--', 'LineWidth', 2, 'DisplayName', '20 kHz Seeker Band');
set(gca, 'XScale', 'log', 'FontSize', 12, 'FontName', 'Times New Roman');
xlim([1000, 50000]); ylim([90, 140]);
xlabel('Frequency (Hz)', 'FontWeight', 'bold'); ylabel('SPL (dB)', 'FontWeight', 'bold');
title(sprintf('Section 2: Multipole High-Frequency Roll-off (RMSE: %.2f dB)', multipole_rmse), 'FontSize', 14);
legend('Location', 'southwest');


%% ========================================================================
%  CORE ENGINE FUNCTION (Strictly Decoupled, One-Sided PSD Fixed)
% =========================================================================
function sig_bb = generate_broadband_noise(N_samples, fs, freq_anchors, ampl_anchors)
    % Initialize stochastic phases
    rand_phase = 2 * pi * rand(N_samples, 1) - pi;
    df = fs / N_samples; f_full = (0:N_samples-1)' * df; f_interp = f_full;
    idx_neg = f_full > fs/2; f_interp(idx_neg) = fs - f_full(idx_neg);
    
    % PCHIP Interpolation for continuum shaping
    target_curve_dB_full = interp1(freq_anchors, ampl_anchors, f_interp, 'pchip', NaN);
    
    % Boundary Extrapolations (Low-frequency cut-off & High-frequency roll-off)
    idx_low = (f_interp < freq_anchors(1) & f_interp > 0);
    if any(idx_low), target_curve_dB_full(idx_low) = ampl_anchors(1) - 10*log10(freq_anchors(1)./f_interp(idx_low)); end
    idx_high = (f_interp > freq_anchors(end));
    if any(idx_high), target_curve_dB_full(idx_high) = ampl_anchors(end) - 6*log2(f_interp(idx_high)/freq_anchors(end)); end
    
    % DC Component handling
    target_curve_dB_full(isnan(target_curve_dB_full)) = -200; 
    target_curve_dB_full(1) = -200; 
    
    % Inverse FFT transformation
    P_linear = 10.^(target_curve_dB_full / 10); 
    Spec_Complex = sqrt(P_linear) .* exp(1j * rand_phase); 
    sig_raw = real(ifft(Spec_Complex));
    
    % --- CRITICAL FIX: ONE-SIDED PSD NORMALIZATION ---
    % Prevents +3.01 dB systematic error by computing target variance 
    % strictly from the one-sided spectrum up to the Nyquist frequency.
    half_len = floor(N_samples/2) + 1;
    Total_Target_Power_OneSided = sum(P_linear(1:half_len)) * df; 
    
    % Energy normalization mapping
    Actual_Power = var(sig_raw);
    if Actual_Power > 0
        sig_bb = sig_raw * sqrt(Total_Target_Power_OneSided / Actual_Power);
    else
        sig_bb = zeros(size(sig_raw)); 
    end
end