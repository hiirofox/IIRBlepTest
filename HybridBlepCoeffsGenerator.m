function result = HybridBlepCoeffsGenerator()
%HYBRIDBLEPCOEFFSGENERATOR Generate fixed-pole IIR+FIR BLEP coefficients.
%
% The exported modal part uses the "modal" lowpass design.  The short FIR
% tables fit the error between a higher-quality "target" design and that
% fixed-pole modal base:
%
%   BLIT  : Hlp * Hhp
%   BLEP  : (Hlp * Hhp - Hhp) / s
%   BLAMP : (Hlp * Hhp - Hhp) / (s * s)
%
% The printed namespace is meant to replace HybridBlep::HybridBlepCoeffs in
% dsp/IIRBlep.h.

    % ---------------- User-facing design parameters ----------------
    modal_order = 8;
    target_order = 24;
    fc_hz = 23500;
    fc_target_hz = 23500;
    Rp = 0.1;
    modal_Rs = 97.5;
    target_Rs = 220.0;
    fs = 48000;
    dc_block_rad = 2.0*pi*20.0;
    dc_block_order = 3;

    fir_size = 16;
    fir_table_size = 24;
    fir_design_mode = 'time-error'; % TableBlep-style continuous 1D error, then polyphase slicing
    fir_signal_len = 65536;
    fir_window_ridge = 1e-4;
    fir_taper = 'right-square-welch'; % 'none', 'right-square-welch', or 'tableblep'
    fir_dc_window = 'square-welch'; % 'square-welch', 'blackman-harris', or 'flat'

    response_samples = 40;
    use_R_tf_export = false;
    minreal_tol = 1e-7;
    % ----------------------------------------------------------------

    cfg = struct();
    cfg.modal_order = modal_order;
    cfg.target_order = target_order;
    cfg.fc_hz = fc_hz;
    cfg.fc_target_hz = fc_target_hz;
    cfg.Rp = Rp;
    cfg.modal_Rs = modal_Rs;
    cfg.target_Rs = target_Rs;
    cfg.fs = fs;
    cfg.dc_block_rad = dc_block_rad;
    cfg.dc_block_order = dc_block_order;
    cfg.fir_size = fir_size;
    cfg.fir_table_size = fir_table_size;
    cfg.fir_design_mode = fir_design_mode;
    cfg.fir_signal_len = fir_signal_len;
    cfg.fir_window_ridge = fir_window_ridge;
    cfg.fir_taper = fir_taper;
    cfg.fir_dc_window = fir_dc_window;
    cfg.response_samples = response_samples;
    cfg.use_R_tf_export = use_R_tf_export;
    cfg.minreal_tol = minreal_tol;
    validate_config(cfg);

    Ts = 1 / fs;
    wc = 2*pi*fc_hz;
    wc_target = 2*pi*fc_target_hz;

    [b_lp_modal, a_lp_modal, Hlp0_modal_before, Hlp0_modal] = ...
        design_elliptic_lowpass(modal_order, Rp, modal_Rs, wc);
    [b_lp_target, a_lp_target, Hlp0_target_before, Hlp0_target] = ...
        design_elliptic_lowpass(target_order, Rp, target_Rs, wc_target);

    [z_hp, p_hp, k_hp] = butter(dc_block_order, dc_block_rad, 'high', 's');
    [b_hp, a_hp] = zp2tf(z_hp, p_hp, k_hp);

    s = tf('s');
    Hhp = tf(b_hp, a_hp);
    Hlp_modal = tf(b_lp_modal, a_lp_modal);
    Hlp_target = tf(b_lp_target, a_lp_target);

    % Readable transfer-function forms.  The polynomial path below is used by
    % default because high-order tf subtraction can leave tiny fake s=0 poles.
    R_modal_blit = minreal(Hlp_modal * Hhp, minreal_tol);
    R_modal_blep = minreal((Hlp_modal * Hhp - Hhp) / s, minreal_tol);
    R_modal_blamp = minreal((Hlp_modal * Hhp - Hhp) / (s * s), minreal_tol);

    R_target_blit = minreal(Hlp_target * Hhp, minreal_tol);
    R_target_blep = minreal((Hlp_target * Hhp - Hhp) / s, minreal_tol);
    R_target_blamp = minreal((Hlp_target * Hhp - Hhp) / (s * s), minreal_tol);

    modal_ba = build_residual_ba( ...
        b_lp_modal, a_lp_modal, b_hp, a_hp, dc_block_order, ...
        {R_modal_blit, R_modal_blep, R_modal_blamp}, use_R_tf_export);
    target_ba = build_residual_ba( ...
        b_lp_target, a_lp_target, b_hp, a_hp, dc_block_order, ...
        {R_target_blit, R_target_blep, R_target_blamp}, use_R_tf_export);

    models(1) = decompose_residual('blit',  modal_ba.blit.b,  modal_ba.blit.a,  1.0 * Ts);
    models(2) = decompose_residual('blep',  modal_ba.blep.b,  modal_ba.blep.a,  1.0);
    models(3) = decompose_residual('blamp', modal_ba.blamp.b, modal_ba.blamp.a, 1.0 / Ts);

    target_models(1) = decompose_residual('blit',  target_ba.blit.b,  target_ba.blit.a,  1.0 * Ts);
    target_models(2) = decompose_residual('blep',  target_ba.blep.b,  target_ba.blep.a,  1.0);
    target_models(3) = decompose_residual('blamp', target_ba.blamp.b, target_ba.blamp.a, 1.0 / Ts);

    [pole_error_blep, pole_error_blamp] = check_poles_are_shared(models);
    fir_tables = design_fir_tables(models, target_models, Ts, cfg);

    fprintf('====================================================\n');
    fprintf('HybridBlep fixed-pole IIR + FIR coefficient export\n');
    fprintf('Modal elliptic order      = %d\n', modal_order);
    fprintf('Target elliptic order     = %d\n', target_order);
    fprintf('Modal cutoff fc           = %.12g Hz\n', fc_hz);
    fprintf('Target cutoff fc          = %.12g Hz\n', fc_target_hz);
    fprintf('Passband ripple Rp        = %.12g dB\n', Rp);
    fprintf('Modal stop Rs             = %.12g dB\n', modal_Rs);
    fprintf('Target stop Rs            = %.12g dB\n', target_Rs);
    fprintf('Sample rate fs            = %.12g Hz\n', fs);
    fprintf('Ts                        = %.12g\n', Ts);
    fprintf('Modal Hlp(0) before/after = %.12g / %.12g\n', real(Hlp0_modal_before), real(Hlp0_modal));
    fprintf('Target Hlp(0) before/after= %.12g / %.12g\n', real(Hlp0_target_before), real(Hlp0_target));
    fprintf('HP cutoff                 = %.12g rad/s\n', dc_block_rad);
    fprintf('HP poles (rad/s)          = ');
    fprintf('%.12g%+.12gj ', [real(p_hp(:)).'; imag(p_hp(:)).']);
    fprintf('\n');
    fprintf('Pole consistency error    = BLEP %.12g, BLAMP %.12g rad/s\n', pole_error_blep, pole_error_blamp);
    fprintf('FIR design                = %s, size %d, table %d, signal %d\n', fir_design_mode, fir_size, fir_table_size, fir_signal_len);
    fprintf('FIR window ridge          = %.12g\n', fir_window_ridge);
    fprintf('FIR one-sided taper       = %s\n', fir_taper);
    fprintf('FIR DC window             = %s\n', fir_dc_window);
    fprintf('Two-pole modal count      = %d\n', numel(models(1).p_two));
    fprintf('One-pole modal count      = %d\n', numel(models(1).p_one));
    for i = 1:numel(models)
        fprintf('%-5s modal h(0+)          = %.12g %+.12gj\n', upper(models(i).name), real(sum(models(i).residues)), imag(sum(models(i).residues)));
        fprintf('%-5s FIR max abs          = %.12g\n', upper(models(i).name), max(abs(fir_tables.(models(i).name)(:))));
    end
    fprintf('====================================================\n\n');

    print_cpp_namespace(Ts, models, fir_tables, cfg);
    plot_hybrid_residuals(models, target_models, fir_tables, fs, cfg);

    result = struct();
    result.cfg = cfg;
    result.fs = fs;
    result.Ts = Ts;
    result.fc_hz = fc_hz;
    result.fc_target_hz = fc_target_hz;
    result.b_lp_modal = b_lp_modal;
    result.a_lp_modal = a_lp_modal;
    result.b_lp_target = b_lp_target;
    result.a_lp_target = a_lp_target;
    result.b_hp = b_hp;
    result.a_hp = a_hp;
    result.z_hp = z_hp;
    result.p_hp = p_hp;
    result.k_hp = k_hp;
    result.models = models;
    result.target_models = target_models;
    result.fir_tables = fir_tables;
    result.pole_error_blep = pole_error_blep;
    result.pole_error_blamp = pole_error_blamp;
end

function validate_config(cfg)
    validateattributes(cfg.modal_order,      {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.target_order,     {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fc_hz,            {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.fc_target_hz,     {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.Rp,               {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.modal_Rs,         {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.target_Rs,        {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.fs,               {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.dc_block_rad,     {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.dc_block_order,   {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fir_size,         {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fir_table_size,   {'numeric'}, {'scalar','integer','>=', 3});
    validateattributes(cfg.fir_signal_len,   {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fir_window_ridge, {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.response_samples, {'numeric'}, {'scalar','integer','positive'});
end

function [b_lp, a_lp, Hlp0_before, Hlp0] = design_elliptic_lowpass(n, Rp, Rs, wc)
    [z, p0, k0] = ellipap(n, Rp, Rs);
    [b0, a0] = zp2tf(z, p0, k0);
    [b_lp, a_lp] = lp2lp(b0, a0, wc);

    Hlp0_before = polyval(b_lp, 0) / polyval(a_lp, 0);
    b_lp = b_lp / Hlp0_before;
    Hlp0 = polyval(b_lp, 0) / polyval(a_lp, 0);
end

function ba = build_residual_ba(b_lp, a_lp, b_hp, a_hp, dc_block_order, R, use_R_tf_export)
    if use_R_tf_export
        [ba.blit.b, ba.blit.a] = tf_to_ba(R{1});
        [ba.blep.b, ba.blep.a] = tf_to_ba(R{2});
        [ba.blamp.b, ba.blamp.a] = tf_to_ba(R{3});
        return;
    end

    b_delta = force_trailing_zeros(poly_sub(b_lp, a_lp), 1, 'Hlp - 1');
    b_hp = force_trailing_zeros(b_hp, dc_block_order, 'Hhp');
    b_shared_residual = force_trailing_zeros(conv(b_delta, b_hp), dc_block_order + 1, '(Hlp - 1) * Hhp');

    ba.blit.b = conv(b_lp, b_hp);
    ba.blit.a = conv(a_lp, a_hp);
    ba.blep.b = b_shared_residual;
    ba.blep.a = conv(conv(a_lp, a_hp), [1 0]);
    ba.blamp.b = b_shared_residual;
    ba.blamp.a = conv(conv(a_lp, a_hp), [1 0 0]);

    [ba.blit.b, ba.blit.a] = cancel_origin_common_factors(ba.blit.b, ba.blit.a);
    [ba.blep.b, ba.blep.a] = cancel_origin_common_factors(ba.blep.b, ba.blep.a);
    [ba.blamp.b, ba.blamp.a] = cancel_origin_common_factors(ba.blamp.b, ba.blamp.a);
end

function fir_tables = design_fir_tables(models, target_models, Ts, cfg)
    fir_tables = struct();
    fir_tables.blit = design_stage_fir_table(models(1), target_models(1), Ts, cfg);
    fir_tables.blep = design_stage_fir_table(models(2), target_models(2), Ts, cfg);
    fir_tables.blamp = design_stage_fir_table(models(3), target_models(3), Ts, cfg);
end

function table = design_stage_fir_table(base_model, target_model, Ts, cfg)
    signal = design_continuous_fir_signal(base_model, target_model, Ts, cfg);
    table = slice_continuous_fir_signal(signal, cfg.fir_size, cfg.fir_table_size);
end

function signal = design_continuous_fir_signal(base_model, target_model, Ts, cfg)
    switch lower(cfg.fir_design_mode)
        case 'time-error'
            sample_offset = linspace(0, cfg.fir_size, cfg.fir_signal_len);
            t = sample_offset * Ts;
            signal = modal_response(target_model, t) - modal_response(base_model, t);
        otherwise
            error('Unknown fir_design_mode "%s". Use "time-error" for the TableBlep-style continuous polyphase table.', cfg.fir_design_mode);
    end

    taper = fir_taper_window(numel(signal), cfg.fir_taper);
    dc_window = dc_compensation_window(numel(signal), cfg.fir_dc_window);
    signal = optimize_windowed_fir_signal(signal, taper, dc_window, 0.0, cfg.fir_window_ridge);
end

function signal = optimize_windowed_fir_signal(desired, taper, dc_window, desired_sum, ridge)
    desired = real(desired(:)).';
    taper = real(taper(:)).';
    dc_window = real(dc_window(:)).';

    if numel(taper) ~= numel(desired) || numel(dc_window) ~= numel(desired)
        error('FIR taper, DC window, and desired signal lengths must match.');
    end

    taper_power = taper .* taper;
    if ridge <= 0
        signal = desired;
        signal(abs(taper) < 1e-12) = 0.0;
        signal = apply_fir_dc_compensation_with_window(signal, desired_sum, dc_window .* (abs(taper) >= 1e-12));
        return;
    end

    fit = taper_power ./ (taper_power + ridge);
    correction_shape = fit .* dc_window;
    correction_sum = sum(correction_shape);

    if abs(correction_sum) < 1e-12
        signal = fit .* desired;
        signal = apply_fir_dc_compensation_with_window(signal, desired_sum, correction_shape);
        return;
    end

    beta = (desired_sum - sum(fit .* desired)) / correction_sum;
    signal = fit .* desired + beta * correction_shape;
end

function table = slice_continuous_fir_signal(signal, fir_size, table_size)
    phase_count = table_size - 1;
    ntable = fir_size * table_size;
    table = zeros(table_size, fir_size);

    for phase = 0:phase_count
        for tap = 0:fir_size-1
            % Match TableBlep: shape one long signal first, then sample it
            % with k = tap*numTables + phase.  The actual terminal zero
            % lives outside the exported polyphase table.
            k = tap * phase_count + phase;
            pos = 1 + (k * (numel(signal) - 1) / ntable);
            table(phase + 1, tap + 1) = interp_vector(signal, pos);
        end
    end
end

function y = interp_vector(x, pos)
    idx = floor(pos);
    frac = pos - idx;
    idx = max(1, min(idx, numel(x) - 1));
    y = x(idx) + frac * (x(idx + 1) - x(idx));
end

function tau = fir_table_tau(index, table_size)
    phase_count = table_size - 1;
    tau = (index - 1) / phase_count;
end

function taper = fir_taper_window(n, taper_name)
    if n <= 1
        taper = ones(1, n);
        return;
    end

    x = (0:n-1) / (n - 1);
    switch lower(taper_name)
        case 'none'
            taper = ones(1, n);
        case 'right-square-welch'
            taper = (1.0 - x.^2).^2;
        case 'tableblep'
            taper = blackman_harris_window(0.5 + 0.5 * x.^10);
        otherwise
            error('Unknown FIR one-sided taper "%s".', taper_name);
    end
end

function c = apply_fir_dc_compensation(c, desired_sum, window_name)
    c = real(c(:)).';
    w = dc_compensation_window(numel(c), window_name);
    c = apply_fir_dc_compensation_with_window(c, desired_sum, w);
end

function c = apply_fir_dc_compensation_with_window(c, desired_sum, w)
    c = real(c(:)).';
    w = real(w(:)).';
    window_sum = sum(w);
    if abs(window_sum) < 1e-12
        w = ones(size(c));
        window_sum = sum(w);
    end

    current_sum = sum(c);
    c = c - (current_sum - desired_sum) * (w / window_sum);
end

function w = dc_compensation_window(n, window_name)
    if n <= 1
        w = ones(1, n);
        return;
    end

    x = (0:n-1) / (n - 1);
    switch lower(window_name)
        case 'square-welch'
            w = (1.0 - (2.0*x - 1.0).^2).^2;
        case 'blackman-harris'
            w = blackman_harris_window(x);
        case 'flat'
            w = ones(1, n);
        otherwise
            error('Unknown FIR DC compensation window "%s".', window_name);
    end
end

function w = blackman_harris_window(x)
    a = 2.0*pi*x;
    w = 0.35875 - 0.48829*cos(a) + 0.14128*cos(2.0*a) - 0.01168*cos(3.0*a);
end

function model = decompose_residual(name, b, a, norm_gain)
    b = trim_leading_zeros(b) * norm_gain;
    a = trim_leading_zeros(a);

    [r_all, p_all, k_dir] = residue(b, a);
    if numel(k_dir) > 1
        warning('%s residual returned a polynomial direct term with %d coefficients. Current C++ modal path only expects scalar direct feedthrough.', upper(name), numel(k_dir));
    end
    if ~isempty(k_dir) && any(abs(k_dir) > 1e-9)
        warning('%s residual has a nonzero direct term. The current C++ HybridBlep modal path does not consume direct gains.', upper(name));
    end

    [p_two, r_two, p_one, r_one] = split_modal_poles(p_all, r_all);

    model = struct();
    model.name = name;
    model.norm_gain = norm_gain;
    model.b = b;
    model.a = a;
    model.poles = p_all;
    model.residues = r_all;
    model.k_dir = k_dir;
    model.directGain = sum(k_dir);
    model.p_two = p_two;
    model.r_two = r_two;
    model.p_one = p_one;
    model.r_one = r_one;
end

function [b, a] = tf_to_ba(H)
    [b, a] = tfdata(H, 'v');
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);
    [b, a] = cancel_origin_common_factors(b, a);
end

function [pole_error_blep, pole_error_blamp] = check_poles_are_shared(models)
    p_ref = modal_pole_vector(models(1));
    p_blep = modal_pole_vector(models(2));
    p_blamp = modal_pole_vector(models(3));

    pole_error_blep = compare_poles(p_ref, p_blep, 'BLIT', 'BLEP');
    pole_error_blamp = compare_poles(p_ref, p_blamp, 'BLIT', 'BLAMP');
end

function err = compare_poles(p_ref, p_test, ref_name, test_name)
    if numel(p_ref) ~= numel(p_test)
        error('%s and %s do not have the same pole count: %d vs %d.', ref_name, test_name, numel(p_ref), numel(p_test));
    end

    if isempty(p_ref)
        err = 0;
        return;
    end

    d = p_ref(:) - p_test(:);
    err = max(abs(d));
    tol = 1e-5 * max(1, max(abs(p_ref(:))));
    if err > tol
        error('%s and %s pole locations differ. max error = %.12g, tol = %.12g.', ref_name, test_name, err, tol);
    end
end

function p = modal_pole_vector(model)
    p = [model.p_two(:); model.p_one(:)];
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

function c = poly_sub(a, b)
    n = max(numel(a), numel(b));
    a = [zeros(1, n - numel(a)), a(:).'];
    b = [zeros(1, n - numel(b)), b(:).'];
    c = a - b;
end

function v = force_trailing_zeros(v, max_count, label)
    v = v(:).';
    scale = max([1, abs(v)]);
    tol = 1e-8 * scale;

    for i = 0:max_count-1
        idx = numel(v) - i;
        if idx < 1
            return;
        end

        if abs(v(idx)) <= tol
            v(idx) = 0;
        else
            warning('%s was expected to have an origin zero at trailing coefficient %d, but abs(coeff)=%.12g exceeds tol=%.12g.', label, i + 1, abs(v(idx)), tol);
            return;
        end
    end
end

function [b, a] = cancel_origin_common_factors(b, a)
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);
    scale = max([1, abs(b), abs(a)]);
    tol = 1e-10 * scale;

    while numel(b) > 1 && numel(a) > 1 && abs(b(end)) <= tol && abs(a(end)) <= tol
        b = b(1:end-1);
        a = a(1:end-1);
    end
end

function v = trim_leading_zeros(v)
    v = v(:).';
    idx = find(v ~= 0, 1, 'first');
    if isempty(idx)
        v = 0;
    else
        v = v(idx:end);
    end
end

function plot_hybrid_residuals(models, target_models, fir_tables, fs, cfg)
    response_samples = cfg.response_samples;
    t_samples = (0:response_samples-1) / fs;
    t_analog = linspace(0, t_samples(end), max(4096, response_samples * 128));

    figure('Name', 'HybridBlep fixed-pole residual impulse responses');
    for i = 1:numel(models)
        fir_table = fir_tables.(models(i).name);
        h_base = modal_response(models(i), t_samples);
        h_target = modal_response(target_models(i), t_samples);
        h_hybrid = h_base;
        copy_count = min(numel(h_hybrid), cfg.fir_size);
        h_hybrid(1:copy_count) = h_hybrid(1:copy_count) + fir_table(1, 1:copy_count);

        h_base_analog = modal_response(models(i), t_analog);
        h_target_analog = modal_response(target_models(i), t_analog);

        subplot(numel(models), 2, 2*i - 1);
        plot(t_analog * 1000, h_target_analog, 'LineWidth', 1.15);
        hold on;
        plot(t_analog * 1000, h_base_analog, '--', 'LineWidth', 0.95);
        plot(t_samples * 1000, h_hybrid, '.', 'MarkerSize', 10);
        plot([t_analog(1), t_analog(end)] * 1000, [0, 0], ':');
        grid on;
        xlabel('Time (ms)');
        ylabel('h[n]');
        title(sprintf('%s tau=0: target analog, base analog, hybrid samples', upper(models(i).name)));
        legend({'target analog', 'base analog', 'hybrid samples'}, 'Location', 'best');
        xlim([0, t_analog(end) * 1000]);

        subplot(numel(models), 2, 2*i);
        plot_flattened_fir_table(fir_table, cfg.fir_table_size, cfg.fir_size);
        grid on;
        xlabel('Sample offset, tap + tau');
        ylabel('FIR coefficient');
        title(sprintf('%s FIR correction table, flattened phase view', upper(models(i).name)));
        xlim([0, cfg.fir_size]);
    end
end

function plot_flattened_fir_table(fir_table, table_size, fir_size)
    phase_count = table_size - 1;
    valid_rows = table_size;
    phase_axis = 0:phase_count;
    hold on;

    x_all = zeros(1, valid_rows * fir_size);
    y_all = zeros(1, valid_rows * fir_size);
    for tap = 0:fir_size-1
        idx = tap * valid_rows + (1:valid_rows);
        k = tap * phase_count + phase_axis;
        x_all(idx) = k / table_size;
        y_all(idx) = fir_table(1:valid_rows, tap + 1).';
    end

    plot(x_all, y_all, 'LineWidth', 1.05);
    plot([0, fir_size], [0, 0], ':', 'LineWidth', 0.8);
    yl = ylim;
    for tap = 1:fir_size-1
        plot([tap, tap], yl, ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.5);
    end
    ylim(yl);
end

function h = modal_response(model, t)
    h = real(sum(bsxfun(@times, model.residues(:), exp(model.poles(:) * t)), 1));
end

function print_cpp_namespace(Ts, models, fir_tables, cfg)
    p_two = models(1).p_two;
    p_one = models(1).p_one;

    fprintf('namespace HybridBlepCoeffs\n');
    fprintf('{\n');
    fprintf('\tconstexpr static float Ts = %.12gf;\n', Ts);
    fprintf('\tconstexpr static int NumTwoPoles = %d;\n', numel(p_two));
    fprintf('\tconstexpr static int NumOnePoles = %d;\n\n', numel(p_one));

    print_cpp_array(two_pole_params_vector(p_two), 'twoPoleParams', 'NumTwoPoles * 2', 2);
    fprintf('\n');
    print_cpp_array(real(p_one(:)).', 'onePoleParams', 'NumOnePoles', 1);
    fprintf('\n');

    print_cpp_array(two_pole_residue_vector(models(1).r_two), 'twoPoleBlitResidues', 'NumTwoPoles * 2', 2);
    print_cpp_array(real(models(1).r_one(:)).', 'onePoleBlitResidues', 'NumOnePoles', 1);
    fprintf('\n');
    print_cpp_array(two_pole_residue_vector(models(2).r_two), 'twoPoleBlepResidues', 'NumTwoPoles * 2', 2);
    print_cpp_array(real(models(2).r_one(:)).', 'onePoleBlepResidues', 'NumOnePoles', 1);
    fprintf('\n');
    print_cpp_array(two_pole_residue_vector(models(3).r_two), 'twoPoleBlampResidues', 'NumTwoPoles * 2', 2);
    print_cpp_array(real(models(3).r_one(:)).', 'onePoleBlampResidues', 'NumOnePoles', 1);
    fprintf('\n');

    fprintf('\tconst float blitDirectGain = %.12g;\n', real(models(1).directGain));
    fprintf('\tconst float blepDirectGain = %.12g;\n', real(models(2).directGain));
    fprintf('\tconst float blampDirectGain = %.12g;\n\n', real(models(3).directGain));

    fprintf('\tconstexpr static int FirSize = %d;\n', cfg.fir_size);
    fprintf('\tconstexpr static int FirTableSize = %d;\n', cfg.fir_table_size);
    print_cpp_matrix(fir_tables.blit, 'firBlitCoeffs', 'FirTableSize', 'FirSize');
    print_cpp_matrix(fir_tables.blep, 'firBlepCoeffs', 'FirTableSize', 'FirSize');
    print_cpp_matrix(fir_tables.blamp, 'firBlampCoeffs', 'FirTableSize', 'FirSize');
    fprintf('}\n');
end

function v = two_pole_params_vector(p_two)
    v = zeros(1, numel(p_two) * 2);
    for i = 1:numel(p_two)
        v((i-1)*2 + 1) = real(p_two(i));
        v((i-1)*2 + 2) = imag(p_two(i));
    end
end

function v = two_pole_residue_vector(r_two)
    v = zeros(1, numel(r_two) * 2);
    for i = 1:numel(r_two)
        v((i-1)*2 + 1) = real(r_two(i));
        v((i-1)*2 + 2) = imag(r_two(i));
    end
end

function print_cpp_array(v, name, size_expr, stride)
    fprintf('\tconst float %s[%s] =\n', name, size_expr);
    fprintf('\t{\n');

    if isempty(v)
        fprintf('\t};\n');
        return;
    end

    for i = 1:stride:numel(v)
        fprintf('\t\t');
        for j = 0:stride-1
            idx = i + j;
            if idx <= numel(v)
                fprintf('%.12gf', v(idx));
            end
            if j + 1 < stride && idx + 1 <= numel(v)
                fprintf(', ');
            end
        end
        if i + stride - 1 < numel(v)
            fprintf(',');
        end
        fprintf('\n');
    end

    fprintf('\t};\n');
end

function print_cpp_matrix(m, name, rows_expr, cols_expr)
    fprintf('\tconst float %s[%s][%s] =\n', name, rows_expr, cols_expr);
    fprintf('\t{\n');
    for r = 1:size(m, 1)
        fprintf('\t\t{ ');
        for c = 1:size(m, 2)
            fprintf('%.12gf', m(r, c));
            if c < size(m, 2)
                fprintf(', ');
            end
        end
        if r < size(m, 1)
            fprintf(' },');
        else
            fprintf(' }');
        end
        fprintf('\n');
    end
    fprintf('\t};\n');
end
