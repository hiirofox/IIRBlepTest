function result = siminphasestepresiduecoeffsgenerator(pole_order, zero_order, fc_hz, fs, fit_oversample, fit_audio_samples, response_samples, fit_method, cep_floor_db, design_pad_factor)
% Minimum-phase Si(x)-1 step-residual modal export.
%
% Idea:
%   1. Reproduce the successful TableBlep target generation:
%      centered Blackman-Harris-windowed sinc -> cepstrum minimum phase ->
%      nonlinear causal tail window -> cumsum / total - 1 -> DC compensation.
%   2. By default, fit residues over a stable fixed analog pole grid. stmcb()
%      and prony() remain available as experimental alternatives.
%   3. Export modal tables directly from the fitted analog poles/residues.
%
% Export format:
%   twoPoleParams = [pre, pim(rad/s), rre, rim, ...]
%   onePoleParams = [pre, rre, ...]
%
% cep_floor_db and design_pad_factor are kept in the signature only for
% compatibility with earlier experiments; the TableBlep-style target uses the
% hard 1e-100 cepstrum floor and cepn = n * 8 from the C++ reference.

    if nargin < 1 || isempty(pole_order),        pole_order = 16; end
    if nargin < 2 || isempty(zero_order),        zero_order = max(0, pole_order - 1); end
    if nargin < 3 || isempty(fc_hz),             fc_hz = 22000; end
    if nargin < 4 || isempty(fs),                fs = 48000; end
    if nargin < 5 || isempty(fit_oversample),    fit_oversample = 10; end
    if nargin < 6 || isempty(fit_audio_samples), fit_audio_samples = 28; end
    if nargin < 7 || isempty(response_samples),  response_samples = 500; end
    if nargin < 8 || isempty(fit_method),        fit_method = 'fixedpoles'; end
    if nargin < 9 || isempty(cep_floor_db),      cep_floor_db = -180; end
    if nargin < 10 || isempty(design_pad_factor), design_pad_factor = 4; end

    validateattributes(pole_order,        {'numeric'}, {'scalar','integer','positive'});
    validateattributes(zero_order,        {'numeric'}, {'scalar','integer','nonnegative'});
    validateattributes(fc_hz,             {'numeric'}, {'scalar','positive'});
    validateattributes(fs,                {'numeric'}, {'scalar','positive'});
    validateattributes(fit_oversample,    {'numeric'}, {'scalar','integer','positive'});
    validateattributes(fit_audio_samples, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(response_samples,  {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cep_floor_db,      {'numeric'}, {'scalar'});
    validateattributes(design_pad_factor, {'numeric'}, {'scalar','integer','positive'});

    if zero_order >= pole_order
        warning('zero_order >= pole_order can create a direct digital term. pole-only modal export will still refit residues without that direct term.');
    end

    if fc_hz >= fs / 2
        warning('fc_hz should stay below the audio Nyquist frequency. Current fc_hz = %.12g Hz, audio Nyquist = %.12g Hz.', fc_hz, fs / 2);
    end

    band_limit = fc_hz / (fs / 2);
    table_window_samples = fit_audio_samples;
    fit_samples = 2^nextpow2(max(65536, table_window_samples * fit_oversample));
    actual_fit_oversample = (fit_samples - 1) / table_window_samples;
    fit_fs = fs * actual_fit_oversample;
    design_samples = fit_samples;
    fit_method = lower(char(fit_method));
    table_cep_multiplier = 8;
    stmcb_iterations = 30;

    table_target = make_tableblep_style_target(fit_samples, table_window_samples, band_limit, fs, table_cep_multiplier);
    t_fit = table_target.t_causal;
    t_design = t_fit;
    t_sinc_centered = table_target.t_centered;
    h_sinc_centered = table_target.sinc_centered;
    h_sinc_windowed = table_target.sinc_windowed;
    h_lp_minphase = table_target.blit_minphase_raw;
    h_lp_area = table_target.blit_sum;
    h_step_minphase = table_target.step_minphase;
    h_residual_design = table_target.blep_after_dc;
    h_residual = h_residual_design;

    if ~all(isfinite(h_residual))
        error('TableBlep-style step residual contains NaN or Inf before IIR fitting.');
    end

    fixed_pole_fit = struct();
    if strcmp(fit_method, 'fixedpoles') || strcmp(fit_method, 'fixed') || strcmp(fit_method, 'fixed-poles')
        [p_all, r_all, fixed_pole_fit] = fit_fixed_analog_modal(h_residual, t_fit, pole_order, fs, fc_hz);
        method_used = 'fixedpoles-ls';
        bd = [];
        ad = [];
        z_poles_raw = exp(p_all / fit_fs);
        z_poles = z_poles_raw;
        stabilized_count = 0;
    else
        [bd, ad, method_used] = fit_iir_impulse(h_residual, zero_order, pole_order, fit_method, stmcb_iterations);
        bd = real(bd(:).');
        ad = real(ad(:).');
        validate_iir_coefficients(bd, ad, method_used);
        bd = bd / ad(1);
        ad = ad / ad(1);
        validate_iir_coefficients(bd, ad, method_used);

        z_poles_raw = roots(ad);
        z_poles = stabilize_z_poles(z_poles_raw);
        stabilized_count = nnz(abs(z_poles_raw) >= 1);

        residues_z = fit_residues_to_samples(z_poles, h_residual);
        p_all = log(z_poles(:)) * fit_fs;
        r_all = residues_z(:);
    end

    finite_mask = isfinite(real(p_all)) & isfinite(imag(p_all)) & isfinite(real(r_all)) & isfinite(imag(r_all));
    p_all = p_all(finite_mask);
    r_all = r_all(finite_mask);
    audio_supernyquist_count = nnz(abs(imag(p_all)) > pi * fs);

    h0_modal = sum(r_all);
    h0_display = real(h0_modal);
    h0_error_from_minus_one = h0_display + 1;
    exportResidueNorm = abs(h0_display);
    if ~isfinite(exportResidueNorm) || exportResidueNorm <= eps
        error('Cannot normalize exported residues; modal h(0+) is %.12g %+.12gj.', real(h0_modal), imag(h0_modal));
    end
    r_all_export = r_all / exportResidueNorm;
    h0_export = sum(r_all_export);
    h0_export_error_from_minus_one = h0_export + 1;

    [p_two, r_two, p_one, r_one] = split_modal_poles(p_all, r_all_export);

    twoPoleParams = zeros(numel(p_two) * 4, 1);
    for i = 1:numel(p_two)
        twoPoleParams((i-1)*4 + 1) = real(p_two(i));
        twoPoleParams((i-1)*4 + 2) = imag(p_two(i));
        twoPoleParams((i-1)*4 + 3) = real(r_two(i));
        twoPoleParams((i-1)*4 + 4) = imag(r_two(i));
    end

    onePoleParams = zeros(numel(p_one) * 2, 1);
    for i = 1:numel(p_one)
        onePoleParams((i-1)*2 + 1) = p_one(i);
        onePoleParams((i-1)*2 + 2) = real(r_one(i));
    end

    h_fit_modal = real(evaluate_modal_time(p_all, r_all_export, t_fit));
    h_target_export = h_residual / exportResidueNorm;
    h_target_design_export = h_residual_design / exportResidueNorm;
    fit_error = h_fit_modal - h_target_export;
    fit_rms_error = sqrt(mean(fit_error.^2));
    fit_peak_error = max(abs(fit_error));

    t_audio = (0:response_samples-1) / fs;
    dense_response_samples = min(numel(t_fit), max(2, round(response_samples * actual_fit_oversample)));
    t_dense_response = t_fit(1:dense_response_samples);
    h_audio = real(evaluate_modal_time(p_all, r_all_export, t_audio));
    h_dense_response = real(evaluate_modal_time(p_all, r_all_export, t_dense_response));
    h_target_audio = interp1(t_fit, h_target_export, t_audio, 'linear', 0);
    h_target_dense_response = h_target_export(1:dense_response_samples);

    f_hz = unique([0, logspace(-3, log10(max(10 * fs / 2, fc_hz * 2)), 16384)]);
    positive_freq_mask = f_hz > 0;
    w = 2*pi*f_hz(positive_freq_mask);
    H_target_residue = evaluate_fir_freq(h_target_export, fit_fs, f_hz);
    H_target_step_view = NaN(size(f_hz));
    H_target_step_view(positive_freq_mask) = H_target_residue(positive_freq_mask) - real(h_target_export(1)) ./ (1i * w);
    H_residue = NaN(size(f_hz));
    H_step_view = NaN(size(f_hz));
    H_residue(positive_freq_mask) = evaluate_modal_freq(p_all, r_all_export, w);
    s_eval = 1i * w;
    H_step_view(positive_freq_mask) = H_residue(positive_freq_mask) - real(h0_export) ./ s_eval;
    H_residue_abs = abs(H_residue);
    H_target_residue_abs = abs(H_target_residue);
    H_residue_peak = max([H_residue_abs(positive_freq_mask), H_target_residue_abs(positive_freq_mask)]);
    H_residue_db_rel = 20*log10(max(H_residue_abs / H_residue_peak, 1e-300));
    H_target_residue_db_rel = 20*log10(max(H_target_residue_abs / H_residue_peak, 1e-300));
    H_step_abs = abs(H_step_view);
    H_target_step_abs = abs(H_target_step_view);
    H_step_peak = max([H_step_abs(positive_freq_mask), H_target_step_abs(positive_freq_mask)]);
    H_step_db_rel = 20*log10(max(H_step_abs / H_step_peak, 1e-300));
    H_target_step_db_rel = 20*log10(max(H_target_step_abs / H_step_peak, 1e-300));

    w_nyquist = 2*pi*(fs/2);
    H_nyquist_residue = evaluate_modal_freq(p_all, r_all_export, w_nyquist);
    H_nyquist_target_residue = evaluate_fir_freq(h_target_export, fit_fs, fs/2);
    H_nyquist_residue_db_rel = 20*log10(max(abs(H_nyquist_residue) / H_residue_peak, 1e-300));
    H_nyquist_target_residue_db_rel = 20*log10(max(abs(H_nyquist_target_residue) / H_residue_peak, 1e-300));
    H_nyquist_step = H_nyquist_residue - real(h0_export) / (1i * w_nyquist);
    H_nyquist_db_rel = 20*log10(max(abs(H_nyquist_step) / H_step_peak, 1e-300));

    fprintf('====================================================\n');
    fprintf('Si minimum-phase step-residual modal export\n');
    fprintf('Prototype                 = centered sinc -> minimum-phase -> integrate -> subtract 1\n');
    fprintf('IIR fit method            = %s\n', method_used);
    fprintf('Pole order                = %d\n', pole_order);
    fprintf('Zero order for fit        = %d\n', zero_order);
    fprintf('Cutoff fc                 = %.12g Hz\n', fc_hz);
    fprintf('Audio fs                  = %.12g Hz\n', fs);
    fprintf('Audio Nyquist             = %.12g Hz\n', fs / 2);
    fprintf('Fit fs                    = %.12g Hz\n', fit_fs);
    fprintf('Requested min oversample  = %d\n', fit_oversample);
    fprintf('Actual fit oversample     = %.12g\n', actual_fit_oversample);
    fprintf('Fit samples               = %d\n', fit_samples);
    fprintf('Table window samples      = %d\n', table_window_samples);
    fprintf('Table bandLimit           = %.12g\n', band_limit);
    fprintf('Cepstrum FFT multiplier   = %d\n', table_cep_multiplier);
    fprintf('Cepstrum abs floor        = 1e-100\n');
    fprintf('Stabilized z poles        = %d\n', stabilized_count);
    fprintf('Poles above audio Nyquist = %d\n', audio_supernyquist_count);
    fprintf('h(0+) before export norm  = %.12g %+.12gj\n', real(h0_modal), imag(h0_modal));
    fprintf('h(0+) - (-1)              = %.12g\n', h0_error_from_minus_one);
    fprintf('C++ residue export norm   = abs(h(0+)) = %.12g\n', exportResidueNorm);
    fprintf('C++ export h(0+)          = %.12g %+.12gj\n', real(h0_export), imag(h0_export));
    fprintf('C++ export h(0+) - (-1)   = %.12g %+.12gj\n', real(h0_export_error_from_minus_one), imag(h0_export_error_from_minus_one));
    fprintf('Fit RMS error             = %.12g\n', fit_rms_error);
    fprintf('Fit peak error            = %.12g\n', fit_peak_error);
    if isfield(fixed_pole_fit, 'relative_residual')
        fprintf('Fixed-pole LS rel error   = %.12g\n', fixed_pole_fit.relative_residual);
        fprintf('Fixed-pole LS cond est    = %.12g\n', fixed_pole_fit.condition_estimate);
    end
    fprintf('Nyquist target raw level  = %.12g dB relative to raw peak\n', H_nyquist_target_residue_db_rel);
    fprintf('Nyquist fit raw level     = %.12g dB relative to raw peak\n', H_nyquist_residue_db_rel);
    fprintf('Nyquist step-view level   = %.12g dB relative to peak\n', H_nyquist_db_rel);
    fprintf('Two-pole modal count      = %d\n', numel(p_two));
    fprintf('One-pole modal count      = %d\n', numel(p_one));
    fprintf('====================================================\n\n');

    print_cpp_vector(twoPoleParams, 'twoPoleParams', 4);
    fprintf('\n');
    print_cpp_vector(onePoleParams, 'onePoleParams', 2);
    fprintf('\nfloat directGain = 0.0f;\n');

    figure('Name', 'Si minimum-phase step residual modal export');

    subplot(4, 1, 1);
    semilogx(f_hz(positive_freq_mask), H_target_residue_db_rel(positive_freq_mask), ':', 'LineWidth', 1.2);
    hold on;
    semilogx(f_hz(positive_freq_mask), H_residue_db_rel(positive_freq_mask), 'LineWidth', 1.2);
    yl = ylim;
    plot([fs/2, fs/2], yl, '--', 'LineWidth', 1.1);
    ylim(yl);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB, peak-normalized)');
    title(sprintf('Raw residual response before/after identification, target Nyquist = %.2f dB', H_nyquist_target_residue_db_rel));
    legend('pre-identification target', 'modal fit', 'Location', 'best');
    xlim([f_hz(find(positive_freq_mask, 1, 'first')), f_hz(end)]);

    subplot(4, 1, 2);
    semilogx(f_hz(positive_freq_mask), H_target_step_db_rel(positive_freq_mask), ':', 'LineWidth', 1.2);
    hold on;
    semilogx(f_hz(positive_freq_mask), H_step_db_rel(positive_freq_mask), 'LineWidth', 1.2);
    yl = ylim;
    plot([fs/2, fs/2], yl, '--', 'LineWidth', 1.1);
    ylim(yl);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB, peak-normalized)');
    title(sprintf('Zero-initial step view: R - h0/s, Nyquist = %.2f dB', H_nyquist_db_rel));
    legend('pre-identification target', 'modal fit', 'Location', 'best');
    xlim([f_hz(find(positive_freq_mask, 1, 'first')), f_hz(end)]);

    subplot(4, 1, 3);
    plot(t_dense_response * 1000, h_target_dense_response, ':', 'LineWidth', 1.2);
    hold on;
    plot(t_dense_response * 1000, h_dense_response, 'LineWidth', 1.2);
    plot(t_audio * 1000, h_target_audio, '.', 'MarkerSize', 5);
    plot(t_audio * 1000, h_audio, '.', 'MarkerSize', 5);
    plot([t_audio(1), t_audio(end)] * 1000, [0, 0], ':');
    plot([t_audio(1), t_audio(end)] * 1000, [-1, -1], '--');
    grid on;
    xlabel('Time (ms)');
    ylabel('h(t)');
    title(sprintf('Oversampled impulse response, first %d audio samples @ %.0f Hz', response_samples, fs));
    legend('minimum-phase target dense', 'modal fit dense', 'target audio samples', 'fit audio samples', 'Location', 'best');
    xlim([0, t_fit(end) * 1000]);

    subplot(4, 1, 4);
    plot(t_dense_response * 1000, h_dense_response - h_target_dense_response, 'LineWidth', 1.2);
    hold on;
    plot([t_audio(1), t_audio(end)] * 1000, [0, 0], ':');
    grid on;
    xlabel('Time (ms)');
    ylabel('fit error');
    title('Modal fit error');
    xlim([0, t_fit(end) * 1000]);

    result = struct();
    result.pole_order = pole_order;
    result.zero_order = zero_order;
    result.fc_hz = fc_hz;
    result.fs = fs;
    result.fit_fs = fit_fs;
    result.fit_oversample = fit_oversample;
    result.actual_fit_oversample = actual_fit_oversample;
    result.fit_samples = fit_samples;
    result.compat_design_pad_factor = design_pad_factor;
    result.design_samples = design_samples;
    result.table_window_samples = table_window_samples;
    result.table_band_limit = band_limit;
    result.table_cep_multiplier = table_cep_multiplier;
    result.fit_method = method_used;
    result.fixed_pole_fit = fixed_pole_fit;
    result.compat_cep_floor_db = cep_floor_db;
    result.bd = bd;
    result.ad = ad;
    result.z_poles_raw = z_poles_raw;
    result.z_poles = z_poles;
    result.stabilized_count = stabilized_count;
    result.audio_supernyquist_count = audio_supernyquist_count;
    result.poles_all_rad = p_all;
    result.residues_all = r_all;
    result.residues_all_export = r_all_export;
    result.two_poles_rad = p_two;
    result.two_residues = r_two;
    result.one_poles_rad = p_one;
    result.one_residues = r_one;
    result.twoPoleParams = twoPoleParams;
    result.onePoleParams = onePoleParams;
    result.h0_modal = h0_modal;
    result.h0_error_from_minus_one = h0_error_from_minus_one;
    result.exportResidueNorm = exportResidueNorm;
    result.h0_export = h0_export;
    result.h0_export_error_from_minus_one = h0_export_error_from_minus_one;
    result.fit_rms_error = fit_rms_error;
    result.fit_peak_error = fit_peak_error;
    result.t_fit = t_fit;
    result.t_design = t_design;
    result.table_target = table_target;
    result.t_sinc_centered = t_sinc_centered;
    result.h_sinc_centered = h_sinc_centered;
    result.h_sinc_windowed = h_sinc_windowed;
    result.h_lp_minphase = h_lp_minphase;
    result.h_lp_area = h_lp_area;
    result.h_step_minphase = h_step_minphase;
    result.h_residual_design = h_residual_design;
    result.h_residual = h_residual;
    result.h_target_export = h_target_export;
    result.h_target_design_export = h_target_design_export;
    result.h_fit_modal = h_fit_modal;
    result.t_audio = t_audio;
    result.t_dense_response = t_dense_response;
    result.h_audio = h_audio;
    result.h_dense_response = h_dense_response;
    result.h_target_audio = h_target_audio;
    result.h_target_dense_response = h_target_dense_response;
    result.f_hz = f_hz;
    result.H_residue = H_residue;
    result.H_target_residue = H_target_residue;
    result.H_residue_db_rel = H_residue_db_rel;
    result.H_target_residue_db_rel = H_target_residue_db_rel;
    result.H_step_view = H_step_view;
    result.H_target_step_view = H_target_step_view;
    result.H_step_db_rel = H_step_db_rel;
    result.H_target_step_db_rel = H_target_step_db_rel;
    result.H_nyquist_step = H_nyquist_step;
    result.H_nyquist_target_residue = H_nyquist_target_residue;
    result.H_nyquist_residue_db_rel = H_nyquist_residue_db_rel;
    result.H_nyquist_target_residue_db_rel = H_nyquist_target_residue_db_rel;
    result.H_nyquist_db_rel = H_nyquist_db_rel;
end

function target = make_tableblep_style_target(n, window_samples, band_limit, fs, cep_multiplier)
    i = 0:n-1;
    normalized = i / n;
    centered = normalized * 2.0 - 1.0;

    input_window = blackman_harris_window(normalized);
    t2 = pi * centered * window_samples / 2.0 * band_limit;
    sinc_values = ones(size(t2));
    nonzero_mask = abs(t2) >= 1e-6;
    sinc_values(nonzero_mask) = sin(t2(nonzero_mask)) ./ t2(nonzero_mask);
    sinc_centered = sinc_values;
    sinc_windowed = input_window .* sinc_values;

    cepn = n * cep_multiplier;
    x2 = zeros(1, cepn);
    x2(1:n) = sinc_windowed;

    X = fft(x2);
    cep = real(ifft(log(abs(X) + 1e-100)));
    cep(2:cepn/2) = 2.0 * cep(2:cepn/2);
    cep(cepn/2 + 2:end) = 0.0;

    H_minphase = exp(fft(cep));
    x2 = real(ifft(H_minphase));
    blit_minphase_raw = x2(1:n);

    causal_normalized = i / (n - 1);
    causal_window = blackman_harris_window(causal_normalized.^10 * 0.5 + 0.5);
    blit_windowed = blit_minphase_raw .* causal_window;
    blit_sum = sum(blit_windowed);
    if ~isfinite(blit_sum) || abs(blit_sum) <= eps
        error('Cannot normalize TableBlep-style BLIT target; sum is %.12g.', blit_sum);
    end

    step_minphase = cumsum(blit_windowed) / blit_sum;
    blep_before_dc = step_minphase - 1.0;
    [blep_after_dc, dc_sum_before, dc_comp_scale, dc_comp_window] = apply_table_dc_compensation(blep_before_dc);

    target = struct();
    target.n = n;
    target.window_samples = window_samples;
    target.band_limit = band_limit;
    target.cep_multiplier = cep_multiplier;
    target.cepn = cepn;
    target.t_centered = centered * window_samples / (2.0 * fs);
    target.t_causal = causal_normalized * window_samples / fs;
    target.input_window = input_window;
    target.sinc_centered = sinc_centered;
    target.sinc_windowed = sinc_windowed;
    target.blit_minphase_raw = blit_minphase_raw;
    target.causal_window = causal_window;
    target.blit_windowed = blit_windowed;
    target.blit_sum = blit_sum;
    target.step_minphase = step_minphase;
    target.blep_before_dc = blep_before_dc;
    target.blep_after_dc = blep_after_dc;
    target.dc_sum_before = dc_sum_before;
    target.dc_comp_scale = dc_comp_scale;
    target.dc_comp_window = dc_comp_window;
end

function [res, currentDCSum, scale, dcWindow] = apply_table_dc_compensation(res)
    n = numel(res);
    currentDCSum = sum(res);
    scale = 0;
    dcWindow = zeros(size(res));

    if n <= 1 || abs(currentDCSum) < 1e-9
        return;
    end

    dcWindow = blackman_harris_window((0:n-1) / (n - 1));
    windowSum = sum(dcWindow);
    if abs(windowSum) < 1e-9
        return;
    end

    scale = currentDCSum / windowSum;
    res = res - dcWindow * scale;
end

function y = blackman_harris_window(x)
    y = zeros(size(x));
    mask = x >= 0 & x <= 1;
    xm = 2.0 * pi * x(mask);
    y(mask) = 0.35875 ...
        - 0.48829 * cos(xm) ...
        + 0.14128 * cos(2.0 * xm) ...
        - 0.01168 * cos(3.0 * xm);
end

function [t, h, h_windowed] = make_windowed_sinc_prototype(n, fs, wc, alpha)
    center = floor(n / 2);
    t = ((0:n-1) - center) / fs;
    h = zeros(size(t));
    nonzero_mask = t ~= 0;
    h(nonzero_mask) = sin(wc * t(nonzero_mask)) ./ (pi * t(nonzero_mask));
    h(~nonzero_mask) = wc / pi;

    h_windowed = h .* make_tukey_window(n, alpha);
end

function w = make_tukey_window(n, alpha)
    if n <= 1
        w = ones(1, n);
        return;
    end

    alpha = max(0, min(1, alpha));
    x = (0:n-1) / (n - 1);

    if alpha <= 0
        w = ones(1, n);
        return;
    end

    if alpha >= 1
        w = 0.5 - 0.5 * cos(2*pi*x);
        return;
    end

    w = ones(1, n);
    left_mask = x < alpha / 2;
    right_mask = x >= 1 - alpha / 2;
    w(left_mask) = 0.5 * (1 + cos(pi * (2*x(left_mask)/alpha - 1)));
    w(right_mask) = 0.5 * (1 + cos(pi * (2*x(right_mask)/alpha - 2/alpha + 1)));
end

function w = make_causal_tail_fade_window(n, fade_fraction)
    w = ones(1, n);
    fade_samples = max(0, min(n - 1, round(n * fade_fraction)));
    if fade_samples <= 1
        return;
    end

    fade = 0.5 * (1 + cos(pi * (0:fade_samples-1) / (fade_samples - 1)));
    w(end-fade_samples+1:end) = fade;
end

function c_min = causal_minphase_lifter(c)
    n = numel(c);
    c = c(:).';
    c_min = zeros(size(c));
    c_min(1) = c(1);

    if mod(n, 2) == 0
        c_min(2:n/2) = 2 * c(2:n/2);
        c_min(n/2 + 1) = c(n/2 + 1);
    else
        c_min(2:(n+1)/2) = 2 * c(2:(n+1)/2);
    end
end

function [bd, ad, method_used] = fit_iir_impulse(h, zero_order, pole_order, fit_method, stmcb_iterations)
    h = h(:);
    if strcmp(fit_method, 'stmcb')
        if exist('stmcb', 'file') == 2
            try
                [bd, ad] = stmcb(h, zero_order, pole_order, stmcb_iterations);
                if iir_coefficients_are_valid(bd, ad)
                    method_used = 'stmcb';
                    return;
                end
                warning('stmcb() returned NaN/Inf or an invalid leading denominator coefficient. Falling back to prony().');
            catch ME
                warning('stmcb() failed: %s. Falling back to prony().', ME.message);
            end
        else
            warning('stmcb() is not available. Falling back to prony().');
        end
    end

    if exist('prony', 'file') == 2
        try
            [bd, ad] = prony(h, zero_order, pole_order);
            if iir_coefficients_are_valid(bd, ad)
                method_used = 'prony';
                return;
            end
            warning('prony() returned NaN/Inf or an invalid leading denominator coefficient. Falling back to LS all-pole fitting.');
        catch ME
            warning('prony() failed: %s. Falling back to LS all-pole fitting.', ME.message);
        end
    end

    [bd, ad] = fit_all_pole_denominator(h, pole_order, zero_order);
    method_used = 'ls-allpole';
end

function [p_all, r_all, info] = fit_fixed_analog_modal(h, t, pole_order, audio_fs, fc_hz)
    h = h(:);
    t = t(:);
    fit_duration = max(t(end) - t(1), 1 / audio_fs);
    p_all = make_fixed_analog_pole_grid(pole_order, audio_fs, fc_hz, fit_duration);

    onset_tau = max(fit_duration * 0.08, 1 / (8 * audio_fs));
    weights = 0.25 + 0.75 * exp(-t / onset_tau);
    weights = weights / max(weights);

    A_time = exp(t * p_all(:).');
    sqrt_weights = sqrt(weights);
    A = bsxfun(@times, A_time, sqrt_weights);
    y = h .* sqrt_weights;

    constraint_weight = sqrt(numel(t)) * 10;
    h0_row = ones(1, numel(p_all));
    area_row = (-1 ./ p_all(:).') / fit_duration;
    A = [A; constraint_weight * h0_row; constraint_weight * area_row];
    y = [y; constraint_weight * h(1); 0];

    col_scale = sqrt(sum(abs(A).^2, 1));
    col_scale(col_scale <= eps) = 1;
    A_scaled = bsxfun(@rdivide, A, col_scale);
    gram = A_scaled' * A_scaled;
    rhs = A_scaled' * y;
    lambda = 1e-10 * max(1, real(trace(gram)) / size(gram, 1));
    r_scaled = (gram + lambda * eye(size(gram))) \ rhs;
    r_all = r_scaled(:) ./ col_scale(:);

    fit = A_time * r_all;
    info = struct();
    info.poles = p_all(:);
    info.weights = weights;
    info.fit_duration = fit_duration;
    info.constraint_weight = constraint_weight;
    info.relative_residual = norm(real(fit) - h) / max(norm(h), eps);
    info.peak_error = max(abs(real(fit) - h));
    info.condition_estimate = cond(gram + lambda * eye(size(gram)));
end

function p_all = make_fixed_analog_pole_grid(pole_order, audio_fs, fc_hz, fit_duration)
    pole_order = max(1, round(pole_order));
    fmax = min(fc_hz * 1.02, audio_fs * 0.499);
    min_decay = 1.25 / fit_duration;
    max_decay = 90.0 / fit_duration;

    real_count = min(pole_order, max(2, round(0.12 * pole_order)));
    if pole_order < 6
        real_count = mod(pole_order, 2);
    end
    if mod(pole_order - real_count, 2) ~= 0
        real_count = real_count + 1;
    end
    real_count = min(real_count, pole_order);
    pair_count = floor((pole_order - real_count) / 2);

    real_poles = [];
    if real_count > 0
        real_poles = -logspace(log10(min_decay), log10(max_decay), real_count);
    end

    complex_poles = [];
    if pair_count > 0
        damping_band_count = min(4, max(1, round(sqrt(pair_count / 2))));
        damping_values = logspace(log10(min_decay * 1.5), log10(max_decay), damping_band_count);
        freq_per_band = ceil(pair_count / damping_band_count);
        freq_values = linspace(fmax / (freq_per_band + 1), fmax, freq_per_band);

        for damping = damping_values
            for freq_hz = freq_values
                if size(complex_poles, 2) >= pair_count
                    break;
                end
                complex_poles(end + 1) = -damping + 1i * 2*pi*freq_hz; %#ok<AGROW>
            end
        end
    end

    p_all = [real_poles(:); complex_poles(:); conj(complex_poles(:))];
    if numel(p_all) > pole_order
        p_all = p_all(1:pole_order);
    end
end

function ok = iir_coefficients_are_valid(bd, ad)
    ok = ~isempty(bd) && ~isempty(ad) && all(isfinite(bd(:))) && all(isfinite(ad(:))) && isfinite(ad(1)) && abs(ad(1)) > eps;
end

function validate_iir_coefficients(bd, ad, method_used)
    if ~iir_coefficients_are_valid(bd, ad)
        error('%s returned invalid IIR coefficients; denominator contains NaN/Inf or ad(1) is zero.', method_used);
    end
end

function [bd, ad] = fit_all_pole_denominator(h, pole_order, zero_order)
    h = h(:);
    first_sample = max(pole_order + 1, zero_order + 2);
    row_count = numel(h) - first_sample + 1;
    if row_count <= pole_order
        error('Not enough samples for LS all-pole fitting.');
    end

    X = zeros(row_count, pole_order);
    y = zeros(row_count, 1);
    for row = 1:row_count
        k = first_sample + row - 1;
        for col = 1:pole_order
            X(row, col) = h(k - col);
        end
        y(row) = -h(k);
    end

    normalizer = max(1, trace(X' * X) / pole_order);
    lambda = 1e-10 * normalizer;
    a_tail = (X' * X + lambda * eye(pole_order)) \ (X' * y);

    ad = [1, a_tail(:).'];
    bd = zeros(1, zero_order + 1);
    bd(1) = h(1);
end

function z = stabilize_z_poles(z)
    z = z(:);
    min_radius = 1e-12;
    max_radius = 0.999999;
    radius = abs(z);

    zero_mask = radius < min_radius;
    z(zero_mask) = min_radius;
    radius = abs(z);

    unstable_mask = radius >= 1;
    z(unstable_mask) = max_radius ./ conj(z(unstable_mask));
end

function residues = fit_residues_to_samples(z_poles, h)
    z_poles = z_poles(:).';
    h = h(:);
    n = (0:numel(h)-1).';
    V = exp(n * log(z_poles));
    residues = V \ h;
end

function h = evaluate_modal_time(poles, residues, t)
    poles = poles(:);
    residues = residues(:);
    t = t(:).';
    h = sum(bsxfun(@times, residues, exp(poles * t)), 1);
end

function H = evaluate_modal_freq(poles, residues, w)
    poles = poles(:);
    residues = residues(:);
    s = 1i * w(:).';
    H = sum(bsxfun(@rdivide, residues, bsxfun(@minus, s, poles)), 1);
end

function H = evaluate_fir_freq(h, sample_rate, f_hz)
    h = h(:).';
    f_hz = f_hz(:).';

    if any(f_hz < 0) || any(f_hz > sample_rate / 2)
        warning('evaluate_fir_freq expects frequencies in [0, sample_rate/2]. Values outside this range will be extrapolated.');
    end

    nfft = 2^nextpow2(max(numel(h) * 2, numel(f_hz) * 4));
    H_fft = fft(h, nfft) / sample_rate;
    positive_count = floor(nfft / 2) + 1;
    f_fft = (0:positive_count-1) * sample_rate / nfft;
    H = interp1(f_fft, H_fft(1:positive_count), f_hz, 'linear', 'extrap');
end

function [p_two, r_two, p_one, r_one] = split_modal_poles(p_all, r_all)
    real_pole_tol = 1e-7;
    real_mask = abs(imag(p_all)) <= real_pole_tol * max(1, abs(p_all));
    pos_imag_mask = imag(p_all) > real_pole_tol * max(1, abs(p_all));

    p_two = p_all(pos_imag_mask);
    r_two = r_all(pos_imag_mask);
    [~, two_order] = sort(imag(p_two));
    p_two = p_two(two_order);
    r_two = r_two(two_order);

    p_one = real(p_all(real_mask));
    r_one = r_all(real_mask);
    [~, one_order] = sort(abs(p_one));
    p_one = p_one(one_order);
    r_one = r_one(one_order);

    if any(abs(imag(r_one)) > 1e-7 * max(1, abs(r_one)))
        warning('Some residues on numerically-real poles are not purely real. Export uses their real parts.');
    end
end

function print_cpp_vector(v, name, stride)
    fprintf('std::vector<float> %s =\n{\n', name);
    if isempty(v)
        fprintf('};\n');
        return;
    end

    for i = 1:stride:numel(v)
        fprintf('    ');
        for j = 0:stride-1
            fprintf('%.12gf', v(i + j));
            if j + 1 < stride
                fprintf(', ');
            end
        end
        if i + stride - 1 < numel(v)
            fprintf(',');
        end
        fprintf('\n');
    end
    fprintf('};\n');
end
