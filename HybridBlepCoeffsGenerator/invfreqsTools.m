function result = HybridBlepCoeffsGenerator3(target_H, n, fs, f_hz)
%HYBRIDBLEPCOEFFSGENERATOR3 Fit complex spectrum samples with analog models.
%
% Usage:
%   result = HybridBlepCoeffsGenerator3()
%   result = HybridBlepCoeffsGenerator3(target_H, n, fs)
%   result = HybridBlepCoeffsGenerator3(target_H, n, fs, f_hz)
%
% target_H is a one-sided complex frequency response sampled at f_hz. If
% f_hz is omitted, samples are assumed to cover 0..fs/2 uniformly. When
% target_H is omitted or empty, a 16384-point stable analog test response is
% generated so the fitting and plots can be exercised immediately.

    cfg = default_config();

    if nargin < 2 || isempty(n),     n = cfg.default_pole_count; end
    if nargin < 3 || isempty(fs),    fs = cfg.default_fs; end
    if nargin < 4,                   f_hz = []; end
    if nargin < 1,                   target_H = []; end

    validateattributes(n,  {'numeric'}, {'scalar','integer','positive'});
    validateattributes(fs, {'numeric'}, {'scalar','positive'});

    test = struct();
    if isempty(target_H)
        [f_hz, target_H, test] = make_test_spectrum(cfg.default_target_points, n, fs);
    else
        target_H = target_H(:);
        if isempty(f_hz)
            f_hz = linspace(0, fs / 2, numel(target_H)).';
        else
            f_hz = f_hz(:);
        end
    end

    validate_spectrum_input(target_H, f_hz);
    w_rad = 2 * pi * f_hz(:);

    invfreqs_model = fit_with_invfreqs(target_H, w_rad, n, cfg);
    tfest_model = fit_with_tfest(target_H, w_rad, n, cfg);
    models = [invfreqs_model, tfest_model];

    fprintf('====================================================\n');
    fprintf('HybridBlepCoeffsGenerator3 analog spectrum fit\n');
    fprintf('Target samples            = %d\n', numel(target_H));
    fprintf('Analog pole count n       = %d\n', n);
    fprintf('Sample rate fs            = %.12g Hz\n', fs);
    fprintf('Target frequency range    = %.12g .. %.12g Hz\n', f_hz(1), f_hz(end));
    fprintf('Analog plot extension     = %d x target frequency range\n', cfg.plot_extend_factor);
    if ~isempty(fieldnames(test))
        fprintf('Target source             = generated stable analog test response\n');
    else
        fprintf('Target source             = user supplied complex spectrum\n');
    end
    fprintf('====================================================\n\n');

    for i = 1:numel(models)
        print_model_summary(models(i));
    end

    plot_data = plot_fit_comparison(f_hz, target_H, models, cfg);

    result = struct();
    result.cfg = cfg;
    result.fs = fs;
    result.n = n;
    result.f_hz = f_hz;
    result.w_rad = w_rad;
    result.target_H = target_H;
    result.test = test;
    result.invfreqs = invfreqs_model;
    result.tfest = tfest_model;
    result.plot = plot_data;
end

function cfg = default_config()
    cfg = struct();
    cfg.default_fs = 48000;
    cfg.default_pole_count = 8;
    cfg.default_target_points = 16384;
    cfg.numerator_order_offset = -1;
    cfg.invfreqs_iterations = 30;
    cfg.invfreqs_tolerance = 1e-8;
    cfg.use_relative_weights = true;
    cfg.relative_weight_floor = 1e-4;
    cfg.plot_extend_factor = 16;
    cfg.plot_line_width = 1.2;
    cfg.pole_real_tol = 1e-8;
end

function validate_spectrum_input(target_H, f_hz)
    validateattributes(target_H, {'numeric'}, {'vector','nonempty'});
    validateattributes(f_hz, {'numeric'}, {'vector','nonempty','real','finite','nonnegative'});

    if numel(target_H) ~= numel(f_hz)
        error('target_H and f_hz must have the same number of elements.');
    end
    if any(~isfinite(real(target_H))) || any(~isfinite(imag(target_H)))
        error('target_H must contain only finite complex values.');
    end
    if any(diff(f_hz(:)) <= 0)
        error('f_hz must be strictly increasing.');
    end
end

function [f_hz, H, test] = make_test_spectrum(num_points, n, fs)
    f_hz = linspace(0, fs / 2, num_points).';

    pair_count = floor(n / 2);
    has_real_pole = mod(n, 2) ~= 0;

    if pair_count > 0
        freq_frac = linspace(0.10, 0.86, pair_count).';
        zeta = linspace(0.08, 0.26, pair_count).';
        wn = 2 * pi * (fs / 2) * freq_frac;
        wd = wn .* sqrt(max(1 - zeta.^2, 0));
        p_upper = -zeta .* wn + 1i * wd;

        residue_phase = linspace(0.35, 1.45, pair_count).';
        residue_scale = (0.75 ./ (1:pair_count).') .* abs(p_upper);
        r_upper = residue_scale .* exp(1i * residue_phase);
    else
        p_upper = zeros(0, 1);
        r_upper = zeros(0, 1);
    end

    if has_real_pole
        p_real = -2 * pi * fs * 0.17;
        r_real = -p_real * 0.35;
    else
        p_real = zeros(0, 1);
        r_real = zeros(0, 1);
    end

    p_all = [p_upper; conj(p_upper); p_real(:)];
    r_all = [r_upper; conj(r_upper); r_real(:)];

    H0 = evaluate_modal_response(r_all, p_all, 0, 0);
    if abs(H0) <= eps
        error('Generated test response has near-zero DC gain.');
    end

    r_all = r_all / H0;
    r_upper = r_upper / H0;
    r_real = r_real / H0;

    H = evaluate_modal_response(r_all, p_all, 0, 2 * pi * f_hz);
    H = H(:);

    test = struct();
    test.poles_all_rad = p_all;
    test.residues_all = r_all;
    test.upper_poles_rad = p_upper;
    test.upper_residues = r_upper;
    test.real_poles_rad = p_real(:);
    test.real_residues = r_real(:);
    test.direct_term = 0;
end

function model = fit_with_invfreqs(H_target, w_rad, n, cfg)
    model = blank_model('invfreqs');

    if ~function_exists('invfreqs')
        model.error_message = 'invfreqs is not available. It requires Signal Processing Toolbox.';
        warning('%s', model.error_message);
        return;
    end

    nb = max(n + cfg.numerator_order_offset, 0);
    na = n;
    weight = make_fit_weights(H_target, cfg);
    w_scale = make_frequency_scale(w_rad);
    w_fit = w_rad(:) / w_scale;

    try
        [b_fit, a_fit] = invfreqs(H_target(:), w_fit(:), nb, na, weight(:), ...
            cfg.invfreqs_iterations, cfg.invfreqs_tolerance);
        [b, a] = denormalize_analog_tf(b_fit, a_fit, w_scale);
        model = package_tf_model(model, b, a, H_target, w_rad, cfg);
        model.frequency_scale_rad_per_sec = w_scale;
    catch ex
        model.error_message = ex.message;
        warning('invfreqs fit failed: %s', ex.message);
    end
end

function model = fit_with_tfest(H_target, w_rad, n, cfg)
    model = blank_model('tfest');

    if ~function_exists('tfest') || ~function_exists('idfrd')
        model.error_message = 'tfest/idfrd is not available. It requires System Identification Toolbox.';
        warning('%s', model.error_message);
        return;
    end

    nz = max(n + cfg.numerator_order_offset, 0);
    w_scale = make_frequency_scale(w_rad);
    w_fit = w_rad(:) / w_scale;

    try
        response_data = reshape(H_target(:), 1, 1, []);
        data = idfrd(response_data, w_fit(:), 0);

        opt = [];
        if function_exists('tfestOptions')
            try
                opt = tfestOptions('Display', 'off');
                if isprop(opt, 'EnforceStability')
                    opt.EnforceStability = true;
                end
            catch
                opt = [];
            end
        end

        if isempty(opt)
            sys = tfest(data, n, nz);
        else
            sys = tfest(data, n, nz, opt);
        end

        [b, a] = tfdata(sys, 'v');
        if iscell(b), b = b{1}; end
        if iscell(a), a = a{1}; end

        [b, a] = denormalize_analog_tf(b, a, w_scale);
        model = package_tf_model(model, b, a, H_target, w_rad, cfg);
        model.sys_normalized = sys;
        if function_exists('tf')
            try
                model.sys = tf(b, a);
            catch
                model.sys = [];
            end
        end
        model.frequency_scale_rad_per_sec = w_scale;
    catch ex
        model.error_message = ex.message;
        warning('tfest fit failed: %s', ex.message);
    end
end

function w_scale = make_frequency_scale(w_rad)
    w_positive = abs(w_rad(w_rad ~= 0));
    if isempty(w_positive)
        w_scale = 1;
    else
        w_scale = max(w_positive);
    end
end

function [b, a] = denormalize_analog_tf(b_fit, a_fit, w_scale)
    b = denormalize_analog_poly(b_fit, w_scale);
    a = denormalize_analog_poly(a_fit, w_scale);
end

function c = denormalize_analog_poly(c_fit, w_scale)
    c_fit = c_fit(:).';
    degree = numel(c_fit) - 1;
    powers = degree:-1:0;
    c = c_fit ./ (w_scale .^ powers);
end

function weight = make_fit_weights(H_target, cfg)
    weight = ones(size(H_target(:)));
    if ~cfg.use_relative_weights
        return;
    end

    peak = max(abs(H_target(:)));
    if peak <= eps
        return;
    end

    floor_mag = peak * cfg.relative_weight_floor;
    weight = 1 ./ max(abs(H_target(:)), floor_mag);
    weight = weight / median(weight);
end

function model = blank_model(name)
    model = struct();
    model.name = name;
    model.valid = false;
    model.error_message = '';
    model.b = [];
    model.a = [];
    model.poles_all_rad = [];
    model.residues_all = [];
    model.direct_terms = [];
    model.upper_poles_rad = [];
    model.upper_residues = [];
    model.real_poles_rad = [];
    model.real_residues = [];
    model.fit_H = [];
    model.rms_abs_error = NaN;
    model.rms_relative_error = NaN;
    model.max_abs_error = NaN;
    model.sys = [];
    model.sys_normalized = [];
    model.frequency_scale_rad_per_sec = NaN;
end

function model = package_tf_model(model, b, a, H_target, w_rad, cfg)
    b = trim_leading_zeros(b(:).');
    a = trim_leading_zeros(a(:).');

    if isempty(a) || a(1) == 0
        error('%s returned an invalid denominator.', model.name);
    end

    b = b / a(1);
    a = a / a(1);

    [r_all, p_all, k_dir] = residue(b, a);
    [p_upper, r_upper, p_real, r_real] = split_upper_modes(p_all, r_all, cfg.pole_real_tol);

    H_fit = evaluate_tf_response(b, a, w_rad);
    err = H_fit(:) - H_target(:);
    peak = max(abs(H_target(:)));
    rel_floor = max(peak * 1e-9, eps);

    model.valid = true;
    model.b = b;
    model.a = a;
    model.poles_all_rad = p_all(:);
    model.residues_all = r_all(:);
    model.direct_terms = k_dir(:);
    model.upper_poles_rad = p_upper(:);
    model.upper_residues = r_upper(:);
    model.real_poles_rad = p_real(:);
    model.real_residues = r_real(:);
    model.fit_H = H_fit(:);
    model.rms_abs_error = sqrt(mean(abs(err).^2));
    model.rms_relative_error = sqrt(mean((abs(err) ./ max(abs(H_target(:)), rel_floor)).^2));
    model.max_abs_error = max(abs(err));
end

function [p_upper, r_upper, p_real, r_real] = split_upper_modes(p_all, r_all, real_tol_factor)
    if isempty(p_all)
        p_upper = [];
        r_upper = [];
        p_real = [];
        r_real = [];
        return;
    end

    tol = real_tol_factor * max(1, max(abs(p_all(:))));
    real_mask = abs(imag(p_all(:))) <= tol;
    upper_mask = imag(p_all(:)) > tol;

    p_upper = p_all(upper_mask);
    r_upper = r_all(upper_mask);
    [~, upper_order] = sort(imag(p_upper));
    p_upper = p_upper(upper_order);
    r_upper = r_upper(upper_order);

    p_real = real(p_all(real_mask));
    r_real = r_all(real_mask);
    [~, real_order] = sort(abs(p_real));
    p_real = p_real(real_order);
    r_real = r_real(real_order);
end

function y = evaluate_tf_response(b, a, w_rad)
    s = 1i * w_rad(:).';
    y = polyval(b, s) ./ polyval(a, s);
    y = y(:);
end

function y = evaluate_modal_response(r, p, k_dir, w_rad)
    s = 1i * w_rad(:).';
    y = zeros(size(s)) + sum(k_dir);
    for i = 1:numel(p)
        y = y + r(i) ./ (s - p(i));
    end
    y = y(:);
end

function plot_data = plot_fit_comparison(f_hz, H_target, models, cfg)
    n_target = numel(f_hz);
    f_ext = linspace(f_hz(1), f_hz(end) * cfg.plot_extend_factor, ...
        max(n_target * cfg.plot_extend_factor, n_target)).';

    H_ext = evaluate_models_on_grid(models, 2 * pi * f_ext);

    figure('Name', 'HybridBlep3 magnitude response');
    hold on;
    plot(f_hz, mag_db(H_target), 'k.', 'MarkerSize', 4);
    legend_entries = {'target'};
    legend_entries = plot_model_family(f_ext, H_ext, models, cfg, legend_entries, 'mag');
    add_target_edge_marker(f_hz(end));
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Magnitude response comparison');
    legend(legend_entries, 'Location', 'best');

    figure('Name', 'HybridBlep3 phase response');
    hold on;
    plot(f_hz, unwrap(angle(H_target)), 'k.', 'MarkerSize', 4);
    legend_entries = {'target'};
    legend_entries = plot_model_family(f_ext, H_ext, models, cfg, legend_entries, 'phase');
    add_target_edge_marker(f_hz(end));
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Unwrapped phase (rad)');
    title('Phase response comparison');
    legend(legend_entries, 'Location', 'best');

    plot_data = struct();
    plot_data.f_extended_hz = f_ext;
    plot_data.H_extended = H_ext;
end

function H_grid = evaluate_models_on_grid(models, w_rad)
    H_grid = cell(size(models));
    for i = 1:numel(models)
        if models(i).valid
            H_grid{i} = evaluate_tf_response(models(i).b, models(i).a, w_rad);
        else
            H_grid{i} = [];
        end
    end
end

function legend_entries = plot_model_family(f_hz, H_grid, models, cfg, legend_entries, mode)
    colors = lines(max(numel(models), 1));
    for i = 1:numel(models)
        if ~models(i).valid
            continue;
        end

        if strcmp(mode, 'mag')
            y = mag_db(H_grid{i});
        else
            y = unwrap(angle(H_grid{i}));
        end

        plot(f_hz, y, 'LineWidth', cfg.plot_line_width, 'Color', colors(i, :));
        legend_entries{end + 1} = models(i).name; %#ok<AGROW>
    end
end

function add_target_edge_marker(f_target_end)
    yl = ylim;
    plot([f_target_end, f_target_end], yl, 'k--', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    ylim(yl);
end

function y = mag_db(H)
    y = 20 * log10(max(abs(H(:)), 1e-300));
end

function print_model_summary(model)
    fprintf('----------------------------------------------------\n');
    fprintf('%s fit\n', model.name);
    if ~model.valid
        fprintf('Status                    = unavailable/failed\n');
        fprintf('Reason                    = %s\n\n', model.error_message);
        return;
    end

    fprintf('Status                    = ok\n');
    fprintf('RMS abs error             = %.12g\n', model.rms_abs_error);
    fprintf('RMS relative error        = %.12g\n', model.rms_relative_error);
    fprintf('Max abs error             = %.12g\n', model.max_abs_error);
    fprintf('Upper conjugate pole count= %d\n', numel(model.upper_poles_rad));
    fprintf('Real pole count           = %d\n', numel(model.real_poles_rad));
    fprintf('Direct term count         = %d\n', numel(model.direct_terms));

    fprintf('\nUpper conjugate poles/residues only:\n');
    if isempty(model.upper_poles_rad)
        fprintf('(none)\n');
    else
        for i = 1:numel(model.upper_poles_rad)
            p = model.upper_poles_rad(i);
            r = model.upper_residues(i);
            fprintf(['  %2d: pole = %.12g %+.12gj rad/s, ', ...
                     'residue = %.12g %+.12gj\n'], ...
                i, real(p), imag(p), real(r), imag(r));
        end
    end

    if ~isempty(model.real_poles_rad)
        fprintf('\nReal poles/residues:\n');
        for i = 1:numel(model.real_poles_rad)
            p = model.real_poles_rad(i);
            r = model.real_residues(i);
            fprintf('  %2d: pole = %.12g rad/s, residue = %.12g %+.12gj\n', ...
                i, p, real(r), imag(r));
        end
    end

    if any(abs(model.direct_terms) > 1e-10)
        fprintf('\nDirect terms from residue():\n');
        for i = 1:numel(model.direct_terms)
            k = model.direct_terms(i);
            fprintf('  %2d: %.12g %+.12gj\n', i, real(k), imag(k));
        end
    end

    fprintf('\n');
end

function y = trim_leading_zeros(x)
    x = x(:).';
    if isempty(x)
        y = x;
        return;
    end

    first = find(x ~= 0, 1, 'first');
    if isempty(first)
        y = 0;
    else
        y = x(first:end);
    end
end

function tf = function_exists(name)
    tf = exist(name, 'file') ~= 0 || exist(name, 'builtin') ~= 0;
end
