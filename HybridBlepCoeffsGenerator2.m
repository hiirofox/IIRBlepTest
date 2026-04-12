function result = HybridBlepCoeffsGenerator2()
%HYBRIDBLEPCOEFFSGENERATOR2 Generate TableBlep FIR + IIR Nyquist correction.
%
% The FIR side intentionally follows dsp/TableBlep.cpp:
%   seed sinc -> minimum phase by real cepstrum -> take first FIR window ->
%   apply TableBlep one-sided window -> BLIT/BLEP DC compensation -> slice.
%
% The target is the same minimum-phase response before the one-sided FIR
% truncation/window.  target_oversamplex controls how much longer the target
% is than the FIR working signal.

    % ---------------- User-facing design parameters ----------------
    fs = 48000;
    band_limit = 0.95;

    fir_size = 16;
    fir_table_size = 24;
    tableblep_signal_len = 65536;
    target_oversamplex = 100;
    use_min_phase = true;

    modal_order = 8;
    fc_hz = band_limit * fs * 0.5;
    Rp = 0.1;
    modal_Rs = 97.5;
    dc_block_rad = 2.0*pi*20.0;
    dc_block_order = 3;

    pole_fc_multipliers = [0.72, 0.82, 0.92, 1.00, 1.08];
    pole_fc_refine = true;
    pole_refine_max_evals = 32;
    freq_fit_start_nyquist = 0.82;
    freq_fit_end_nyquist = 1.0;
    freq_fit_points = 4096;
    freq_fit_ridge = 1e-8;
    freq_fit_dc_weight = 10.0;
    freq_fit_weight_power = 2.0;

    fir_target_plot_points = fir_size * fir_table_size;
    freq_plot_points = 2048;
    % ----------------------------------------------------------------

    cfg = struct();
    cfg.fs = fs;
    cfg.band_limit = band_limit;
    cfg.fir_size = fir_size;
    cfg.fir_table_size = fir_table_size;
    cfg.tableblep_signal_len = tableblep_signal_len;
    cfg.target_oversamplex = target_oversamplex;
    cfg.use_min_phase = use_min_phase;
    cfg.modal_order = modal_order;
    cfg.fc_hz = fc_hz;
    cfg.Rp = Rp;
    cfg.modal_Rs = modal_Rs;
    cfg.dc_block_rad = dc_block_rad;
    cfg.dc_block_order = dc_block_order;
    cfg.pole_fc_multipliers = pole_fc_multipliers;
    cfg.pole_fc_refine = pole_fc_refine;
    cfg.pole_refine_max_evals = pole_refine_max_evals;
    cfg.freq_fit_start_nyquist = freq_fit_start_nyquist;
    cfg.freq_fit_end_nyquist = freq_fit_end_nyquist;
    cfg.freq_fit_points = freq_fit_points;
    cfg.freq_fit_ridge = freq_fit_ridge;
    cfg.freq_fit_dc_weight = freq_fit_dc_weight;
    cfg.freq_fit_weight_power = freq_fit_weight_power;
    cfg.fir_target_plot_points = fir_target_plot_points;
    cfg.freq_plot_points = freq_plot_points;
    validate_config(cfg);

    Ts = 1 / fs;
    data = make_tableblep_fir_and_target(cfg);
    [pole_set, pole_search] = select_shared_poles_from_blep(data.blep, Ts, cfg);

    models(1) = fit_stage_iir_frequency(data.blit,  pole_set, Ts, cfg);
    models(2) = fit_stage_iir_frequency(data.blep,  pole_set, Ts, cfg);
    models(3) = fit_stage_iir_frequency(data.blamp, pole_set, Ts, cfg);

    fir_tables = struct();
    fir_tables.blit = data.blit.fir_table;
    fir_tables.blep = data.blep.fir_table;
    fir_tables.blamp = data.blamp.fir_table;

    fprintf('====================================================\n');
    fprintf('HybridBlep2 TableBlep FIR + fixed-pole IIR frequency correction export\n');
    fprintf('Sample rate fs            = %.12g Hz\n', fs);
    fprintf('Ts                        = %.12g\n', Ts);
    fprintf('Band limit                = %.12g * Nyquist\n', band_limit);
    fprintf('FIR size/table            = %d / %d\n', fir_size, fir_table_size);
    fprintf('TableBlep signal len      = %d\n', tableblep_signal_len);
    fprintf('Target oversamplex        = %d\n', target_oversamplex);
    fprintf('Target high-res len       = %d\n', tableblep_signal_len * target_oversamplex);
    fprintf('Min phase                 = %d\n', use_min_phase);
    fprintf('Modal elliptic order      = %d\n', modal_order);
    fprintf('Modal cutoff fc           = %.12g Hz\n', fc_hz);
    fprintf('Passband ripple Rp        = %.12g dB\n', Rp);
    fprintf('Modal stop Rs             = %.12g dB\n', modal_Rs);
    fprintf('HP cutoff                 = %.12g rad/s\n', dc_block_rad);
    fprintf('Pole fc multipliers       = ');
    fprintf('%.12g ', pole_fc_multipliers);
    fprintf('\n');
    fprintf('Pole fc refine            = %d, max evals = %d\n', pole_fc_refine, pole_refine_max_evals);
    fprintf('Pole fc candidates(score) = ');
    for i = 1:numel(pole_search.candidate_fc_hz)
        fprintf('%.6g:%.6g ', pole_search.candidate_fc_hz(i), pole_search.candidate_scores(i));
    end
    fprintf('\n');
    if pole_search.used_refinement
        fprintf('Pole fc refined(score)    = %.12g:%.12g\n', pole_search.refined_fc_hz, pole_search.refined_score);
    end
    fprintf('Selected BLEP pole fc     = %.12g Hz\n', pole_search.best_fc_hz);
    fprintf('BLEP pole-search score    = %.12g\n', pole_search.best_score);
    fprintf('Freq fit Nyquist range    = %.12g..%.12g, %d points\n', ...
        freq_fit_start_nyquist, freq_fit_end_nyquist, freq_fit_points);
    fprintf('Freq fit ridge/DC/power   = %.12g / %.12g / %.12g\n', ...
        freq_fit_ridge, freq_fit_dc_weight, freq_fit_weight_power);
    fprintf('Two-pole modal count      = %d\n', numel(pole_set.p_two));
    fprintf('One-pole modal count      = %d\n', numel(pole_set.p_one));
    for i = 1:numel(models)
        stage = data.(models(i).name);
        name = upper(models(i).name);
        fprintf('%-5s target tau0 sum      = %.12g\n', name, stage.target_tau0_sum);
        fprintf('%-5s FIR tau0 sum         = %.12g\n', name, stage.fir_tau0_sum);
        fprintf('%-5s IIR sample sum       = %.12g\n', name, iir_sample_sum(models(i), Ts));
        fprintf('%-5s freq err before/after= %.12g / %.12g\n', name, models(i).freq_error_before_rms, models(i).freq_error_after_rms);
        fprintf('%-5s freq err improvement = %.12g dB\n', name, models(i).freq_improvement_db);
        fprintf('%-5s Nyq abs before/after = %.12g / %.12g\n', name, models(i).nyquist_abs_before, models(i).nyquist_abs_after);
        fprintf('%-5s FIR table max abs    = %.12g\n', name, max(abs(stage.fir_table(:))));
    end
    fprintf('====================================================\n\n');

    print_cpp_namespace(Ts, models, fir_tables, cfg);
    plot_hybrid2_debug(models, data, Ts, cfg);

    result = struct();
    result.cfg = cfg;
    result.fs = fs;
    result.Ts = Ts;
    result.data = data;
    result.pole_set = pole_set;
    result.pole_search = pole_search;
    result.models = models;
    result.fir_tables = fir_tables;
end

function validate_config(cfg)
    validateattributes(cfg.fs,                     {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.band_limit,             {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.fir_size,               {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fir_table_size,         {'numeric'}, {'scalar','integer','>=', 3});
    validateattributes(cfg.tableblep_signal_len,   {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.target_oversamplex,     {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.modal_order,            {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.fc_hz,                  {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.Rp,                     {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.modal_Rs,               {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.dc_block_rad,           {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.dc_block_order,         {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.pole_fc_multipliers,    {'numeric'}, {'vector','positive'});
    validateattributes(cfg.pole_refine_max_evals,  {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.freq_fit_start_nyquist, {'numeric'}, {'scalar','>=',0,'<=',1});
    validateattributes(cfg.freq_fit_end_nyquist,   {'numeric'}, {'scalar','>=',0,'<=',1});
    validateattributes(cfg.freq_fit_points,        {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.freq_fit_ridge,         {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.freq_fit_dc_weight,     {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.freq_fit_weight_power,  {'numeric'}, {'scalar','nonnegative'});
    validateattributes(cfg.fir_target_plot_points, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.freq_plot_points,       {'numeric'}, {'scalar','integer','positive'});

    if cfg.freq_fit_start_nyquist >= cfg.freq_fit_end_nyquist
        error('freq_fit_start_nyquist must be smaller than freq_fit_end_nyquist.');
    end
    if ~(islogical(cfg.pole_fc_refine) || (isnumeric(cfg.pole_fc_refine) && isscalar(cfg.pole_fc_refine)))
        error('pole_fc_refine must be a scalar logical or numeric flag.');
    end
end

function data = make_tableblep_fir_and_target(cfg)
    n = cfg.tableblep_signal_len;
    target_n = n * cfg.target_oversamplex;

    % This block mirrors TableBlep::ApplySincBlep() before minimum-phase.
    i = 0:n-1;
    x = i / n;
    t = x * 2.0 - 1.0;
    wd = blackman_harris_window(x);
    arg = pi * t * cfg.fir_size * 0.5 * cfg.band_limit;
    sinc_like = ones(1, n);
    mask = abs(arg) >= 1e-12;
    sinc_like(mask) = sin(arg(mask)) ./ arg(mask);
    seed = wd .* sinc_like;

    % TableBlep used cepn = n*8 and copied the first n samples.  Here the
    % full cepstral buffer is kept as the target, and the FIR still truncates
    % the first n samples exactly like TableBlep does.
    if cfg.use_min_phase
        target_blit_seed = minimum_phase_from_padded_seed(seed, target_n);
    else
        target_blit_seed = zeros(1, target_n);
        target_blit_seed(1:n) = seed;
    end

    causal_taper = blackman_harris_window(x.^10 * 0.5 + 0.5);
    fir_blit_seed = target_blit_seed(1:n) .* causal_taper;

    target_int = sum(target_blit_seed);
    fir_int = sum(fir_blit_seed);
    if abs(target_int) < 1e-18 || abs(fir_int) < 1e-18
        error('TableBlep sinc normalization failed because the integral is too small.');
    end

    phase_norm = cfg.fir_table_size; % TableBlep: numTables + 1.
    target_blit = target_blit_seed / target_int * phase_norm;
    target_blep = cumsum(target_blit_seed) / target_int - 1.0;

    fir_blit = fir_blit_seed / fir_int * phase_norm;
    fir_blep = cumsum(fir_blit_seed) / fir_int - 1.0;

    % Exactly like TableBlep: BLIT/BLEP get DC compensation after shaping;
    % BLAMP is generated from the compensated BLEP and is not compensated.
    fir_blit = apply_dc_compensation(fir_blit, 0.0, 'blackman-harris');
    fir_blep = apply_dc_compensation(fir_blep, 0.0, 'blackman-harris');

    dtv = cfg.fir_size / n;
    target_blamp = cumsum(target_blep) * dtv;
    fir_blamp = cumsum(fir_blep) * dtv;

    data = struct();
    data.offset_target = (0:target_n-1) / n * cfg.fir_size;
    data.offset_fir = (0:n-1) / n * cfg.fir_size;
    data.blit = make_stage_data('blit', target_blit, fir_blit, true, cfg);
    data.blep = make_stage_data('blep', target_blep, fir_blep, false, cfg);
    data.blamp = make_stage_data('blamp', target_blamp, fir_blamp, false, cfg);
end

function stage = make_stage_data(name, target_signal, fir_signal, has_impulse_delta, cfg)
    fir_table = slice_tableblep_signal(fir_signal, cfg.fir_size, cfg.fir_table_size);
    if has_impulse_delta
        fir_table(:, 1) = fir_table(:, 1) - 1.0;
    end

    stage = struct();
    stage.name = name;
    stage.target_signal = target_signal;
    stage.fir_signal = fir_signal;
    stage.fir_table = fir_table;
    stage.has_impulse_delta = has_impulse_delta;
    stage.target_tau0_sum = estimate_target_tau0_sum(target_signal, has_impulse_delta, cfg);
    stage.fir_tau0_sum = sum(fir_table(1, :));
end

function y = minimum_phase_from_padded_seed(seed, fft_len)
    x2 = zeros(1, fft_len);
    x2(1:numel(seed)) = seed;

    X = fft(x2);
    cep = ifft(log(abs(X) + 1e-100));

    % 1-based version of TableBlep's causal cepstrum window:
    % cep[0] keep, cep[1..N/2-1] double, cep[N/2] keep, upper half zero.
    cep(2:fft_len/2) = cep(2:fft_len/2) * 2.0;
    cep(fft_len/2 + 2:end) = 0.0;

    Y = exp(fft(cep));
    y = real(ifft(Y));
end

function [best_pole_set, pole_search] = select_shared_poles_from_blep(blep_stage, Ts, cfg)
    candidate_fc = unique(max(1.0, min(cfg.fs * 0.499, cfg.fc_hz * cfg.pole_fc_multipliers(:).')));
    best_score = Inf;
    best_pole_set = [];
    best_model = [];
    candidate_scores = zeros(size(candidate_fc));
    best_index = 1;

    for i = 1:numel(candidate_fc)
        [score, pole_set, model] = score_blep_pole_fc(candidate_fc(i), blep_stage, Ts, cfg);
        candidate_scores(i) = score;
        if score < best_score
            best_score = score;
            best_pole_set = pole_set;
            best_model = model;
            best_index = i;
        end
    end

    refined_fc = NaN;
    refined_score = NaN;
    used_refinement = false;
    if cfg.pole_fc_refine && numel(candidate_fc) >= 2
        lo = candidate_fc(max(1, best_index - 1));
        hi = candidate_fc(min(numel(candidate_fc), best_index + 1));
        if hi > lo
            opts = optimset('Display', 'off', ...
                'TolX', max(1e-3, best_pole_set.fc_hz * 1e-6), ...
                'MaxFunEvals', cfg.pole_refine_max_evals);
            objective = @(fc) score_blep_pole_fc(fc, blep_stage, Ts, cfg);
            try
                [refined_fc, refined_score] = fminbnd(objective, lo, hi, opts);
                if refined_score < best_score
                    [best_score, best_pole_set, best_model] = score_blep_pole_fc(refined_fc, blep_stage, Ts, cfg);
                    used_refinement = true;
                end
            catch ME
                warning('BLEP pole fc refinement failed: %s', ME.message);
            end
        end
    end

    pole_search = struct();
    pole_search.candidate_fc_hz = candidate_fc;
    pole_search.candidate_scores = candidate_scores;
    pole_search.refined_fc_hz = refined_fc;
    pole_search.refined_score = refined_score;
    pole_search.used_refinement = used_refinement;
    pole_search.best_fc_hz = best_pole_set.fc_hz;
    pole_search.best_score = best_score;
    pole_search.best_blep_model = best_model;
end

function [score, pole_set, model] = score_blep_pole_fc(fc_hz, blep_stage, Ts, cfg)
    pole_set = design_shared_iir_poles_for_fc(cfg, fc_hz);
    model = fit_stage_iir_frequency(blep_stage, pole_set, Ts, cfg);
    score = model.freq_error_after_rms;
end

function pole_set = design_shared_iir_poles_for_fc(cfg, fc_hz)
    wc = 2*pi*fc_hz;
    [z_lp, p_lp, k_lp] = ellipap(cfg.modal_order, cfg.Rp, cfg.modal_Rs);
    [b0, a0] = zp2tf(z_lp, p_lp, k_lp);
    [b_lp, a_lp] = lp2lp(b0, a0, wc);
    Hlp0 = polyval(b_lp, 0) / polyval(a_lp, 0);
    b_lp = b_lp / Hlp0;

    [z_hp, p_hp, k_hp] = butter(cfg.dc_block_order, cfg.dc_block_rad, 'high', 's');
    [b_hp, a_hp] = zp2tf(z_hp, p_hp, k_hp);

    p_all = roots(conv(a_lp, a_hp));
    [p_two, p_one] = split_poles(p_all);

    pole_set = struct();
    pole_set.fc_hz = fc_hz;
    pole_set.p_two = p_two;
    pole_set.p_one = p_one;
    pole_set.p_all = [p_two(:); conj(p_two(:)); p_one(:)];
    pole_set.b_lp = b_lp;
    pole_set.a_lp = a_lp;
    pole_set.z_hp = z_hp;
    pole_set.p_hp = p_hp;
    pole_set.k_hp = k_hp;
    pole_set.b_hp = b_hp;
    pole_set.a_hp = a_hp;
end

function [p_two, p_one] = split_poles(p_all)
    tol = 1e-7;
    real_mask = abs(imag(p_all)) <= tol * max(1, abs(p_all));
    two_mask = imag(p_all) > tol * max(1, abs(p_all));

    p_two = p_all(two_mask);
    [~, two_order] = sort(imag(p_two));
    p_two = p_two(two_order);

    p_one = real(p_all(real_mask));
    [~, one_order] = sort(abs(p_one));
    p_one = p_one(one_order);
end

function model = fit_stage_iir_frequency(stage, pole_set, Ts, cfg)
    w = linspace(cfg.freq_fit_start_nyquist, cfg.freq_fit_end_nyquist, cfg.freq_fit_points).' * pi;
    H_target = target_frequency_response(stage, cfg, w);
    H_fir = fir_frequency_response(stage, w);
    H_want = H_target - H_fir;

    A_complex = modal_frequency_design_matrix(pole_set.p_two, pole_set.p_one, Ts, w);
    weight = frequency_fit_weights(w, cfg);
    Aw_complex = bsxfun(@times, A_complex, weight(:));
    yw_complex = H_want(:) .* weight(:);

    A = [real(Aw_complex); imag(Aw_complex)];
    y = [real(yw_complex); imag(yw_complex)];

    if cfg.freq_fit_dc_weight > 0
        dc_row = modal_sample_sum_row(pole_set.p_two, pole_set.p_one, Ts);
        desired_iir_sum = stage.target_tau0_sum - stage.fir_tau0_sum;
        A = [A; cfg.freq_fit_dc_weight * dc_row];
        y = [y; cfg.freq_fit_dc_weight * desired_iir_sum];
    end

    lhs = A.' * A + cfg.freq_fit_ridge * eye(size(A, 2));
    rhs = A.' * y;
    x = lhs \ rhs;

    [r_two, r_one] = unpack_modal_solution(x, pole_set.p_two, pole_set.p_one);
    p_all = [pole_set.p_two(:); conj(pole_set.p_two(:)); pole_set.p_one(:)];
    r_all = [r_two(:); conj(r_two(:)); r_one(:)];

    H_iir = A_complex * x;
    H_hybrid = H_fir + H_iir;

    model = struct();
    model.name = stage.name;
    model.p_two = pole_set.p_two;
    model.p_one = pole_set.p_one;
    model.r_two = r_two;
    model.r_one = r_one;
    model.poles = p_all;
    model.residues = r_all;
    model.directGain = 0;
    model.has_impulse_delta = stage.has_impulse_delta;
    model.freq_error_before_rms = sqrt(mean(abs(H_fir - H_target).^2));
    model.freq_error_after_rms = sqrt(mean(abs(H_hybrid - H_target).^2));
    model.freq_improvement_db = 20.0 * log10(max(model.freq_error_before_rms, 1e-18) / max(model.freq_error_after_rms, 1e-18));
    model.nyquist_abs_before = abs(H_fir(end));
    model.nyquist_abs_after = abs(H_hybrid(end));
end

function A = modal_frequency_design_matrix(p_two, p_one, Ts, w)
    w = w(:);
    E = exp(-1i * w);
    cols = numel(p_two) * 2 + numel(p_one);
    A = zeros(numel(w), cols);

    col = 1;
    for i = 1:numel(p_two)
        lambda = exp(p_two(i) * Ts);
        q1 = 1.0 ./ (1.0 - lambda * E);
        q2 = 1.0 ./ (1.0 - conj(lambda) * E);
        A(:, col) = q1 + q2;
        A(:, col + 1) = 1i * (q1 - q2);
        col = col + 2;
    end

    for i = 1:numel(p_one)
        lambda = exp(p_one(i) * Ts);
        A(:, col) = 1.0 ./ (1.0 - lambda * E);
        col = col + 1;
    end
end

function row = modal_sample_sum_row(p_two, p_one, Ts)
    row = zeros(1, numel(p_two) * 2 + numel(p_one));

    col = 1;
    for i = 1:numel(p_two)
        q = 1.0 / (1.0 - exp(p_two(i) * Ts));
        row(col) = 2.0 * real(q);
        row(col + 1) = -2.0 * imag(q);
        col = col + 2;
    end

    for i = 1:numel(p_one)
        q = 1.0 / (1.0 - exp(p_one(i) * Ts));
        row(col) = real(q);
        col = col + 1;
    end
end

function [r_two, r_one] = unpack_modal_solution(x, p_two, p_one)
    r_two = zeros(numel(p_two), 1);
    r_one = zeros(numel(p_one), 1);

    col = 1;
    for i = 1:numel(p_two)
        r_two(i) = x(col) + 1i * x(col + 1);
        col = col + 2;
    end

    for i = 1:numel(p_one)
        r_one(i) = x(col);
        col = col + 1;
    end
end

function weight = frequency_fit_weights(w, cfg)
    x = (w(:) / pi - cfg.freq_fit_start_nyquist) / max(1e-12, cfg.freq_fit_end_nyquist - cfg.freq_fit_start_nyquist);
    x = max(0.0, min(1.0, x));
    weight = 0.1 + 0.9 * x .^ cfg.freq_fit_weight_power;
end

function H = target_frequency_response(stage, cfg, w)
    h = target_tau0_sequence(stage, cfg);
    H = freqz(h(:), 1, w(:));
end

function H = fir_frequency_response(stage, w)
    H = freqz(stage.fir_table(1, :).', 1, w(:));
end

function h = target_tau0_sequence(stage, cfg)
    max_sample = floor((numel(stage.target_signal) - 1) / cfg.tableblep_signal_len * cfg.fir_size);
    offsets = 0:max_sample;
    h = target_residual_at(stage, offsets, cfg);
end

function y = target_at(signal, offset_samples, cfg)
    pos = 1 + offset_samples / cfg.fir_size * cfg.tableblep_signal_len;
    y = interp_vector_or_zero(signal, pos);
end

function y = fir_at(stage, offset_samples, cfg)
    pos = 1 + offset_samples / cfg.fir_size * cfg.tableblep_signal_len;
    y = interp_vector_or_zero(stage.fir_signal, pos);
    if stage.has_impulse_delta
        zero_mask = abs(offset_samples) < 1e-12;
        y(zero_mask) = y(zero_mask) - 1.0;
    end
end

function y = target_residual_at(stage, offset_samples, cfg)
    y = target_at(stage.target_signal, offset_samples, cfg);
    if stage.has_impulse_delta
        zero_mask = abs(offset_samples) < 1e-12;
        y(zero_mask) = y(zero_mask) - 1.0;
    end
end

function y = interp_vector_or_zero(x, pos)
    y = zeros(size(pos));
    mask = pos >= 1 & pos <= numel(x);
    if ~any(mask)
        return;
    end

    p = pos(mask);
    idx = floor(p);
    frac = p - idx;

    end_mask = idx >= numel(x);
    idx(end_mask) = numel(x) - 1;
    frac(end_mask) = 1.0;

    idx = max(1, min(idx, numel(x) - 1));
    y(mask) = x(idx) + frac .* (x(idx + 1) - x(idx));
end

function table = slice_tableblep_signal(signal, fir_size, table_size)
    num_tables = table_size - 1;
    ntable = fir_size * table_size;
    table = zeros(table_size, fir_size);

    for phase = 0:num_tables
        for tap = 0:fir_size-1
            k = tap * num_tables + phase;
            pos = 1 + (k * (numel(signal) - 1) / ntable);
            table(phase + 1, tap + 1) = interp_vector_or_zero(signal, pos);
        end
    end
end

function total = estimate_target_tau0_sum(signal, has_impulse_delta, cfg)
    max_sample = floor((numel(signal) - 1) / cfg.tableblep_signal_len * cfg.fir_size);
    offsets = 0:max_sample;
    total = sum(target_at(signal, offsets, cfg));
    if has_impulse_delta
        total = total - 1.0;
    end
end

function s = iir_sample_sum(model, Ts)
    s = 0.0;
    for i = 1:numel(model.p_two)
        lambda = exp(model.p_two(i) * Ts);
        s = s + 2.0 * real(model.r_two(i) / (1.0 - lambda));
    end
    for i = 1:numel(model.p_one)
        lambda = exp(model.p_one(i) * Ts);
        s = s + real(model.r_one(i) / (1.0 - lambda));
    end
end

function h = modal_response(model, t)
    h = modal_response_from_parts(model.p_two, model.r_two, model.p_one, model.r_one, t);
end

function H = modal_frequency_response(model, w, Ts)
    E = exp(-1i * w(:));
    H = zeros(size(E));

    for i = 1:numel(model.p_two)
        lambda = exp(model.p_two(i) * Ts);
        r = model.r_two(i);
        H = H + r ./ (1.0 - lambda * E) + conj(r) ./ (1.0 - conj(lambda) * E);
    end

    for i = 1:numel(model.p_one)
        lambda = exp(model.p_one(i) * Ts);
        H = H + model.r_one(i) ./ (1.0 - lambda * E);
    end
end

function y = amp_db(x)
    y = 20.0 * log10(max(abs(x), 1e-14));
end

function h = modal_response_from_parts(p_two, r_two, p_one, r_one, t)
    t = t(:).';
    h = zeros(size(t));

    for i = 1:numel(p_two)
        h = h + 2.0 * real(r_two(i) * exp(p_two(i) * t));
    end
    for i = 1:numel(p_one)
        h = h + real(r_one(i) * exp(p_one(i) * t));
    end
end

function c = apply_dc_compensation(c, desired_sum, window_name)
    c = real(c(:)).';
    w = dc_compensation_window(numel(c), window_name);
    w = real(w(:)).';
    window_sum = sum(w);
    if abs(window_sum) < 1e-12
        return;
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
            error('Unknown DC compensation window "%s".', window_name);
    end
end

function w = blackman_harris_window(x)
    a = 2.0*pi*x;
    w = 0.35875 - 0.48829*cos(a) + 0.14128*cos(2.0*a) - 0.01168*cos(3.0*a);
end

function plot_hybrid2_debug(models, data, Ts, cfg)
    stages = {data.blit, data.blep, data.blamp};
    fir_target_offset = (0:cfg.fir_target_plot_points-1) / cfg.fir_table_size;
    freq_nyq = linspace(cfg.freq_fit_start_nyquist, cfg.freq_fit_end_nyquist, cfg.freq_plot_points);
    w = freq_nyq(:) * pi;

    figure('Name', 'HybridBlep2 TableBlep FIR and IIR Nyquist correction');
    for i = 1:numel(models)
        stage = stages{i};

        fir_line = fir_at(stage, fir_target_offset, cfg);
        target_line = target_residual_at(stage, fir_target_offset, cfg);

        H_target = target_frequency_response(stage, cfg, w);
        H_fir = fir_frequency_response(stage, w);
        H_iir = modal_frequency_response(models(i), w, Ts);
        H_hybrid = H_fir + H_iir;

        subplot(numel(models), 2, 2*i - 1);
        plot(fir_target_offset, target_line, 'LineWidth', 1.15);
        hold on;
        plot(fir_target_offset, fir_line, '--', 'LineWidth', 1.05);
        plot([fir_target_offset(1), fir_target_offset(end)], [0, 0], ':');
        grid on;
        xlabel('Sample offset');
        ylabel('h');
        title(sprintf('%s target and TableBlep FIR, %d points', upper(stage.name), cfg.fir_target_plot_points));
        legend({'target', 'TableBlep FIR'}, 'Location', 'best');
        xlim([fir_target_offset(1), fir_target_offset(end)]);

        subplot(numel(models), 2, 2*i);
        plot(freq_nyq, amp_db(H_target), 'LineWidth', 1.15);
        hold on;
        plot(freq_nyq, amp_db(H_hybrid), '-', 'LineWidth', 1.05);
        plot(freq_nyq, amp_db(H_fir), '-.', 'LineWidth', 0.95);
        plot(freq_nyq, amp_db(H_iir), '--', 'LineWidth', 1.05);
        grid on;
        xlabel('Frequency / Nyquist');
        ylabel('Magnitude (dB)');
        title(sprintf('%s Nyquist-band frequency correction', upper(stage.name)));
        legend({'target', 'hybrid', 'TableBlep FIR', 'IIR correction'}, 'Location', 'best');
        xlim([freq_nyq(1), freq_nyq(end)]);
    end
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

    fprintf('\tconst float blitDirectGain = 0;\n');
    fprintf('\tconst float blepDirectGain = 0;\n');
    fprintf('\tconst float blampDirectGain = 0;\n\n');

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
