function result = DesignHybridSystem2(n_poles, cfgOverride)
%DESIGNHYBRIDSYSTEM2 IIR-first hybrid BLIT/BLEP/BLAMP design.
%
% This script intentionally uses the parallel hybrid structure:
%
%   target_long[n] ~= iir[n] + fir_correction[n]
%
% The IIR starts at sample 0 and keeps running. The short FIR correction also
% starts at sample 0, but is truncated to firWindowSize * firTableNum samples.

    if nargin < 1 || isempty(n_poles),      n_poles = 8; end
    if nargin < 2 || isempty(cfgOverride),  cfgOverride = struct(); end

    validateattributes(n_poles, {'numeric'}, {'scalar','integer','positive'});

    cfg = default_config();
    cfg = apply_overrides(cfg, cfgOverride);
    cfg = finalize_config(cfg);

    print_config(cfg, n_poles);

    source = call_without_figures(@() MinphaseBandlimitSiGenerator( ...
        cfg.n_target, cfg.bandlimit, cfg.minphaseOversample));
    target = make_target_residuals(source, cfg);
    target_spectra = make_target_spectra(target, cfg);

    fits = fit_iir_to_targets(target_spectra, n_poles, cfg);
    iir = sample_iir_models(fits, target, cfg);

    fir = make_fir_corrections(target, iir, cfg);
    hybrid = make_hybrid_result(target, iir, fir, cfg);

    plots = struct();
    plots.iir_magnitude = plot_iir_magnitude(target_spectra, iir, cfg);
    plots.hybrid_magnitude = plot_hybrid_magnitude(target_spectra, hybrid, cfg);
    plots.hybrid_time = plot_hybrid_time(target, iir, fir, hybrid, cfg);

    result = struct();
    result.cfg = cfg;
    result.source = source;
    result.target = target;
    result.target_spectra = target_spectra;
    result.fits = fits;
    result.iir = iir;
    result.fir = fir;
    result.hybrid = hybrid;
    result.plots = plots;
    result.poles_residues = collect_pole_residue_summary(fits, target.names);
end

function cfg = default_config()
    cfg = struct();
    cfg.fs = 48000;
    cfg.Ts = [];

    cfg.firTableNum = 64;
    cfg.targetSampleRate = [];
    cfg.n_target = 131072 * 8;

    cfg.firWindowSize = 16;
    cfg.firTotalLength = [];

    cfg.cutoffSafety = 0.95;
    cfg.targetCutoffHz = [];
    cfg.bandlimit = [];

    cfg.minphaseOversample = 8;
    cfg.causalWindowPower = 10.0;
    cfg.blitScale = [];
    cfg.blampDtv = [];

    cfg.maxFitPoints = 16384;
    cfg.relativeFitWeightFloor = 1e-4;
    cfg.hideInvfreqsToolFigures = true;
    cfg.runInvfreqsFit = true;
    cfg.enforceStablePoles = true;
    cfg.enforceIirH0Zero = true;
    cfg.minPoleDampingRad = [];

    cfg.iirFitBandwidthMultiplier = 8;
    cfg.iirFitMaxHz = [];

    cfg.applyFirCorrectionWindow = true;
    cfg.applyFirDcCompensation = true;

    cfg.spectrumPlotMaxHz = [];
    cfg.spectrumPlotPoints = 4096;
    cfg.spectrumMagnitudeFloor = 1e-300;
    cfg.timePlotSamples = [];
end

function cfg = apply_overrides(cfg, override)
    validateattributes(override, {'struct'}, {'scalar'});
    names = fieldnames(override);
    for i = 1:numel(names)
        cfg.(names{i}) = override.(names{i});
    end
end

function cfg = finalize_config(cfg)
    validateattributes(cfg.fs, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.firTableNum, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.n_target, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.firWindowSize, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.cutoffSafety, {'numeric'}, {'scalar','positive','<',1});
    validateattributes(cfg.minphaseOversample, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.causalWindowPower, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.maxFitPoints, {'numeric'}, {'scalar','integer','>=',2});
    validateattributes(cfg.relativeFitWeightFloor, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.iirFitBandwidthMultiplier, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.spectrumPlotPoints, {'numeric'}, {'scalar','integer','>=',2});
    validateattributes(cfg.spectrumMagnitudeFloor, {'numeric'}, {'scalar','positive'});

    if mod(cfg.n_target, 2) ~= 0
        error('n_target must be even so the one-sided FFT spectrum has an exact Nyquist bin.');
    end

    cfg.Ts = 1 / cfg.fs;
    cfg.firTotalLength = cfg.firWindowSize * cfg.firTableNum;
    if cfg.firTotalLength > cfg.n_target
        error('firTotalLength (%d) must not exceed n_target (%d).', cfg.firTotalLength, cfg.n_target);
    end

    if isempty(cfg.targetSampleRate)
        cfg.targetSampleRate = cfg.fs * cfg.firTableNum;
    end
    validateattributes(cfg.targetSampleRate, {'numeric'}, {'scalar','positive'});

    if isempty(cfg.targetCutoffHz)
        cfg.targetCutoffHz = (cfg.fs / 2) * cfg.cutoffSafety;
    end
    validateattributes(cfg.targetCutoffHz, {'numeric'}, {'scalar','positive'});
    if cfg.targetCutoffHz >= cfg.targetSampleRate / 2
        error('targetCutoffHz must be below the target-sample-rate Nyquist frequency.');
    end

    if isempty(cfg.bandlimit)
        cfg.bandlimit = 2 * cfg.targetCutoffHz / cfg.targetSampleRate;
    end
    validateattributes(cfg.bandlimit, {'numeric'}, {'scalar','positive','<',1});

    if isempty(cfg.blitScale)
        cfg.blitScale = cfg.firTableNum;
    end
    validateattributes(cfg.blitScale, {'numeric'}, {'scalar','positive'});

    if isempty(cfg.blampDtv)
        cfg.blampDtv = 1 / cfg.firTableNum;
    end
    validateattributes(cfg.blampDtv, {'numeric'}, {'scalar','positive'});

    if isempty(cfg.minPoleDampingRad)
        cfg.minPoleDampingRad = cfg.targetSampleRate * 1e-6;
    end
    validateattributes(cfg.minPoleDampingRad, {'numeric'}, {'scalar','positive'});

    if isempty(cfg.iirFitMaxHz)
        cfg.iirFitMaxHz = min(cfg.targetSampleRate / 2, ...
            cfg.targetCutoffHz * cfg.iirFitBandwidthMultiplier);
    end
    validateattributes(cfg.iirFitMaxHz, {'numeric'}, {'scalar','positive'});
    cfg.iirFitMaxHz = min(cfg.iirFitMaxHz, cfg.targetSampleRate / 2);

    if isempty(cfg.spectrumPlotMaxHz)
        cfg.spectrumPlotMaxHz = cfg.iirFitMaxHz;
    end
    validateattributes(cfg.spectrumPlotMaxHz, {'numeric'}, {'scalar','positive'});
    cfg.spectrumPlotMaxHz = min(cfg.spectrumPlotMaxHz, cfg.targetSampleRate / 2);

    if isempty(cfg.timePlotSamples)
        cfg.timePlotSamples = min(cfg.n_target, max(cfg.firTotalLength, 2048));
    end
    validateattributes(cfg.timePlotSamples, {'numeric'}, {'scalar','integer','positive'});
    cfg.timePlotSamples = min(cfg.timePlotSamples, cfg.n_target);

    cfg.hideInvfreqsToolFigures = scalar_flag(cfg.hideInvfreqsToolFigures, 'hideInvfreqsToolFigures');
    cfg.runInvfreqsFit = scalar_flag(cfg.runInvfreqsFit, 'runInvfreqsFit');
    cfg.enforceStablePoles = scalar_flag(cfg.enforceStablePoles, 'enforceStablePoles');
    cfg.enforceIirH0Zero = scalar_flag(cfg.enforceIirH0Zero, 'enforceIirH0Zero');
    cfg.applyFirCorrectionWindow = scalar_flag(cfg.applyFirCorrectionWindow, 'applyFirCorrectionWindow');
    cfg.applyFirDcCompensation = scalar_flag(cfg.applyFirDcCompensation, 'applyFirDcCompensation');
end

function y = scalar_flag(x, name)
    if ~(islogical(x) || (isnumeric(x) && isscalar(x)))
        error('%s must be a scalar logical/numeric flag.', name);
    end
    y = logical(x);
end

function print_config(cfg, n_poles)
    fprintf('====================================================\n');
    fprintf('Hybrid BLIT/BLEP/BLAMP system design, IIR first\n');
    fprintf('base fs / Ts              = %.12g Hz / %.12g s\n', cfg.fs, cfg.Ts);
    fprintf('target sample rate        = %.12g Hz\n', cfg.targetSampleRate);
    fprintf('n_target                  = %d\n', cfg.n_target);
    fprintf('target cutoff             = %.12g Hz\n', cfg.targetCutoffHz);
    fprintf('normalized bandlimit      = %.12g\n', cfg.bandlimit);
    fprintf('sinc first zero           = %.12g target samples\n', 1 / cfg.bandlimit);
    fprintf('FIR window/table/total    = %d / %d / %d\n', ...
        cfg.firWindowSize, cfg.firTableNum, cfg.firTotalLength);
    fprintf('IIR pole count            = %d\n', n_poles);
    fprintf('IIR fit max frequency     = %.12g Hz\n', cfg.iirFitMaxHz);
    fprintf('IIR h(0+) constrained zero= %d\n', cfg.enforceIirH0Zero);
    fprintf('FIR correction window     = %d\n', cfg.applyFirCorrectionWindow);
    fprintf('FIR correction DC comp    = %d\n', cfg.applyFirDcCompensation);
    fprintf('BLIT scale / BLAMP dtv    = %.12g / %.12g\n', cfg.blitScale, cfg.blampDtv);
    fprintf('====================================================\n\n');
end

function target = make_target_residuals(source, cfg)
    sinc_unit_area = source.minphase_sinc(:);

    target = struct();
    target.names = {'blit', 'blep', 'blamp'};
    target.blit = sinc_unit_area * cfg.blitScale;
    target.blep = cumsum(sinc_unit_area) - 1.0;

    % Integrate from the sample boundary so BLAMP starts at 0.
    target.blamp = exclusive_cumsum(target.blep) * cfg.blampDtv;
end

function spectra = make_target_spectra(target, cfg)
    one_sided_count = cfg.n_target / 2 + 1;
    f_hz = (0:one_sided_count-1).' * (cfg.targetSampleRate / cfg.n_target);

    fit_all = find(f_hz <= cfg.iirFitMaxHz);
    if isempty(fit_all)
        error('No frequency samples fall inside iirFitMaxHz.');
    end
    local_fit_idx = fit_indices(numel(fit_all), min(cfg.maxFitPoints, numel(fit_all)));
    fit_idx = fit_all(local_fit_idx);

    spectra = struct();
    spectra.names = target.names;
    spectra.freq_hz = f_hz;
    spectra.fit_idx = fit_idx;
    spectra.fit_freq_hz = f_hz(fit_idx);
    spectra.full_fft = struct();
    spectra.fit_H = struct();

    for i = 1:numel(target.names)
        name = target.names{i};
        X = fft(target.(name)(:));
        H = X(1:one_sided_count);
        spectra.full_fft.(name) = H;
        spectra.fit_H.(name) = H(fit_idx);
    end
end

function fits = fit_iir_to_targets(spectra, n_poles, cfg)
    fits = struct();

    if ~cfg.runInvfreqsFit
        fprintf('Skipping invfreqsTools fitting because cfg.runInvfreqsFit is false.\n');
        for i = 1:numel(spectra.names)
            name = spectra.names{i};
            fits.(name) = blank_fit();
        end
        return;
    end

    for i = 1:numel(spectra.names)
        name = spectra.names{i};
        fprintf('Fitting %s target spectrum with an IIR model, 0..%.12g Hz...\n', ...
            upper(name), cfg.iirFitMaxHz);

        H = spectra.fit_H.(name);
        f_hz = spectra.fit_freq_hz;
        tool_fit = call_invfreqs_tools(H, n_poles, cfg.targetSampleRate, f_hz, cfg);

        constrained = struct();
        if isfield(tool_fit, 'invfreqs') && tool_fit.invfreqs.valid
            constrained = refit_residues(tool_fit.invfreqs, H, f_hz, cfg);
        else
            warning('%s invfreqsTools fit did not return a valid invfreqs model.', upper(name));
        end

        fits.(name) = struct();
        fits.(name).tool = tool_fit;
        fits.(name).refit = constrained;
    end
end

function fit = blank_fit()
    fit = struct();
    fit.tool = struct();
    fit.refit = struct();
    fit.refit.valid = false;
end

function fit = call_invfreqs_tools(H, n_poles, sample_rate, f_hz, cfg)
    if cfg.hideInvfreqsToolFigures
        fit = call_without_figures(@() call_invfreqs_tools_visible(H, n_poles, sample_rate, f_hz));
    else
        fit = call_invfreqs_tools_visible(H, n_poles, sample_rate, f_hz);
    end
end

function fit = call_invfreqs_tools_visible(H, n_poles, sample_rate, f_hz)
    try
        fit = invfreqsTools(H, n_poles, sample_rate, f_hz);
    catch ex
        if exist('HybridBlepCoeffsGenerator3', 'file') ~= 0
            fit = HybridBlepCoeffsGenerator3(H, n_poles, sample_rate, f_hz);
        else
            rethrow(ex);
        end
    end
end

function constrained = refit_residues(model, H, f_hz, cfg)
    p_upper = model.upper_poles_rad(:);
    p_real = model.real_poles_rad(:);

    if isempty(p_upper) && isempty(p_real)
        error('Cannot refit residues because the model has no poles.');
    end

    if cfg.enforceStablePoles
        [p_upper, p_real] = force_stable_poles(p_upper, p_real, cfg);
    end

    A = sampled_modal_basis_matrix(f_hz(:), p_upper, p_real, cfg.targetSampleRate, cfg.n_target);
    weights = make_relative_weights(H, cfg.relativeFitWeightFloor);
    Aw = A .* weights(:);
    bw = H(:) .* weights(:);

    Ar = [real(Aw); imag(Aw)];
    br = [real(bw); imag(bw)];

    C = h0_constraint_row(numel(p_upper), numel(p_real));
    if cfg.enforceIirH0Zero
        Z = null(C);
        if isempty(Z)
            theta = zeros(size(C, 2), 1);
        else
            theta = Z * ((Ar * Z) \ br);
        end
    else
        theta = Ar \ br;
    end

    [r_upper, r_real, p_all, r_all] = unpack_modal_params(theta, p_upper, p_real);
    H_fit = A * theta;
    err = H_fit - H(:);
    h0 = C * theta;

    constrained = struct();
    constrained.valid = true;
    constrained.poles_all_rad = p_all;
    constrained.residues_all = r_all;
    constrained.upper_poles_rad = p_upper;
    constrained.upper_residues = r_upper;
    constrained.real_poles_rad = p_real;
    constrained.real_residues = r_real;
    constrained.direct_terms = [];
    constrained.h0 = h0;
    constrained.h0_constrained_zero = cfg.enforceIirH0Zero;
    constrained.fit_domain = 'sampled_modal_finite_dft';
    constrained.sample_rate = cfg.targetSampleRate;
    constrained.sample_count = cfg.n_target;
    constrained.fit_H = H_fit;
    constrained.rms_abs_error = sqrt(mean(abs(err).^2));
    constrained.max_abs_error = max(abs(err));
end

function [p_upper, p_real] = force_stable_poles(p_upper, p_real, cfg)
    p_upper = force_stable_pole_vector(p_upper, cfg.minPoleDampingRad);
    p_real = real(force_stable_pole_vector(p_real, cfg.minPoleDampingRad));
end

function p = force_stable_pole_vector(p, minDampingRad)
    if isempty(p)
        return;
    end

    re = real(p);
    im = imag(p);
    bad = re >= -minDampingRad;
    re(bad) = -max(abs(re(bad)), minDampingRad);
    p = re + 1i * im;
end

function iir = sample_iir_models(fits, target, cfg)
    iir = struct();
    iir.names = target.names;
    iir.time = struct();
    iir.full_fft = struct();
    iir.valid = struct();

    one_sided_count = cfg.n_target / 2 + 1;

    for i = 1:numel(target.names)
        name = target.names{i};
        if isfield(fits, name) && isfield(fits.(name), 'refit') && ...
                isfield(fits.(name).refit, 'valid') && fits.(name).refit.valid
            x = sampled_modal_time_response(fits.(name).refit, cfg.n_target, cfg.targetSampleRate);
            X = fft(x);
            iir.time.(name) = x;
            iir.full_fft.(name) = X(1:one_sided_count);
            iir.valid.(name) = true;
        else
            iir.time.(name) = zeros(cfg.n_target, 1);
            iir.full_fft.(name) = zeros(one_sided_count, 1);
            iir.valid.(name) = false;
        end
    end
end

function x = sampled_modal_time_response(model, n, sampleRate)
    t = (0:n-1).' / sampleRate;
    x = zeros(n, 1);

    p = model.poles_all_rad(:);
    r = model.residues_all(:);
    for i = 1:numel(p)
        x = x + r(i) * exp(p(i) * t);
    end

    if isfield(model, 'direct_terms') && ~isempty(model.direct_terms)
        x(1) = x(1) + sum(model.direct_terms);
    end
    x = real(x);
end

function fir = make_fir_corrections(target, iir, cfg)
    fir = struct();
    fir.names = target.names;
    fir.full_residual = struct();
    fir.raw = struct();
    fir.windowed = struct();
    fir.table_saved = struct();
    fir.padded = struct();
    fir.window = one_sided_blackman_window(cfg.firTotalLength, cfg.causalWindowPower);

    for i = 1:numel(target.names)
        name = target.names{i};

        residual = target.(name)(:) - iir.time.(name)(:);
        raw = residual(1:cfg.firTotalLength);

        if cfg.applyFirCorrectionWindow
            windowed = raw .* fir.window;
        else
            windowed = raw;
        end

        if cfg.applyFirDcCompensation
            saved = apply_dc_compensation(windowed);
        else
            saved = windowed;
        end

        padded = zeros(cfg.n_target, 1);
        padded(1:cfg.firTotalLength) = saved;

        fir.full_residual.(name) = residual;
        fir.raw.(name) = raw;
        fir.windowed.(name) = windowed;
        fir.table_saved.(name) = saved;
        fir.padded.(name) = padded;
    end
end

function hybrid = make_hybrid_result(target, iir, fir, cfg)
    one_sided_count = cfg.n_target / 2 + 1;

    hybrid = struct();
    hybrid.names = target.names;
    hybrid.time = struct();
    hybrid.full_fft = struct();
    hybrid.error_time = struct();
    hybrid.error_fft = struct();
    hybrid.rms_time_error = struct();
    hybrid.max_time_error = struct();

    for i = 1:numel(target.names)
        name = target.names{i};
        x = iir.time.(name) + fir.padded.(name);
        err = x - target.(name)(:);
        X = fft(x);
        E = fft(err);

        hybrid.time.(name) = x;
        hybrid.full_fft.(name) = X(1:one_sided_count);
        hybrid.error_time.(name) = err;
        hybrid.error_fft.(name) = E(1:one_sided_count);
        hybrid.rms_time_error.(name) = sqrt(mean(err .^ 2));
        hybrid.max_time_error.(name) = max(abs(err));
    end
end

function plot_data = plot_iir_magnitude(spectra, iir, cfg)
    f_hz = spectra.freq_hz;
    plot_idx = plot_indices_for_frequency(f_hz, cfg.spectrumPlotMaxHz, cfg.spectrumPlotPoints);

    figure('Name', 'DesignHybridSystem2 IIR-only magnitude fit');
    for i = 1:numel(spectra.names)
        name = spectra.names{i};
        subplot(numel(spectra.names), 1, i);
        plot(f_hz(plot_idx), mag_db(spectra.full_fft.(name)(plot_idx), cfg), 'LineWidth', 1.1);
        hold on;
        plot(f_hz(plot_idx), mag_db(iir.full_fft.(name)(plot_idx), cfg), '--', 'LineWidth', 1.1);
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('%s target vs IIR-only magnitude', upper(name)));
        legend({'target long FIR', 'IIR only'}, 'Location', 'best');
        xlim([f_hz(plot_idx(1)), f_hz(plot_idx(end))]);
    end

    plot_data = struct();
    plot_data.plot_idx = plot_idx;
end

function plot_data = plot_hybrid_magnitude(spectra, hybrid, cfg)
    f_hz = spectra.freq_hz;
    plot_idx = plot_indices_for_frequency(f_hz, cfg.spectrumPlotMaxHz, cfg.spectrumPlotPoints);

    figure('Name', 'DesignHybridSystem2 hybrid magnitude after FIR correction');
    for i = 1:numel(spectra.names)
        name = spectra.names{i};
        subplot(numel(spectra.names), 1, i);
        plot(f_hz(plot_idx), mag_db(spectra.full_fft.(name)(plot_idx), cfg), 'LineWidth', 1.1);
        hold on;
        plot(f_hz(plot_idx), mag_db(hybrid.full_fft.(name)(plot_idx), cfg), '--', 'LineWidth', 1.1);
        plot(f_hz(plot_idx), mag_db(hybrid.error_fft.(name)(plot_idx), cfg), ':', 'LineWidth', 1.0);
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('%s target vs IIR + short FIR correction', upper(name)));
        legend({'target long FIR', 'IIR + short FIR', 'hybrid error'}, 'Location', 'best');
        xlim([f_hz(plot_idx(1)), f_hz(plot_idx(end))]);
    end

    plot_data = struct();
    plot_data.plot_idx = plot_idx;
end

function plot_data = plot_hybrid_time(target, iir, fir, hybrid, cfg)
    idx = (1:cfg.timePlotSamples).';
    x = (idx - 1) / cfg.firTableNum;

    figure('Name', 'DesignHybridSystem2 hybrid time response');
    for i = 1:numel(target.names)
        name = target.names{i};
        subplot(numel(target.names), 1, i);
        plot(x, target.(name)(idx), 'LineWidth', 1.1);
        hold on;
        plot(x, hybrid.time.(name)(idx), '--', 'LineWidth', 1.1);
        plot(x, iir.time.(name)(idx), ':', 'LineWidth', 1.0);
        plot(x, fir.padded.(name)(idx), '-.', 'LineWidth', 1.0);
        grid on;
        xlabel('Base-rate sample offset');
        ylabel(name);
        title(sprintf('%s target, hybrid, IIR, short FIR correction', upper(name)));
        legend({'target long FIR', 'IIR + short FIR', 'IIR only', 'short FIR correction'}, ...
            'Location', 'best');
    end

    plot_data = struct();
    plot_data.sample_offset = x;
end

function summary = collect_pole_residue_summary(fits, names)
    summary = struct();
    for i = 1:numel(names)
        name = names{i};
        item = struct();
        if isfield(fits, name) && isfield(fits.(name), 'refit') && ...
                isfield(fits.(name).refit, 'valid') && fits.(name).refit.valid
            model = fits.(name).refit;
            item.valid = true;
            item.upper_poles_rad = model.upper_poles_rad;
            item.upper_residues = model.upper_residues;
            item.real_poles_rad = model.real_poles_rad;
            item.real_residues = model.real_residues;
            item.h0 = model.h0;
            item.rms_abs_error = model.rms_abs_error;
            item.max_abs_error = model.max_abs_error;
        else
            item.valid = false;
            item.upper_poles_rad = [];
            item.upper_residues = [];
            item.real_poles_rad = [];
            item.real_residues = [];
            item.h0 = [];
            item.rms_abs_error = [];
            item.max_abs_error = [];
        end
        summary.(name) = item;
    end
end

function A = sampled_modal_basis_matrix(f_hz, p_upper, p_real, sampleRate, sampleCount)
    pair_count = numel(p_upper);
    real_count = numel(p_real);
    A = zeros(numel(f_hz), 2 * pair_count + real_count);

    for i = 1:pair_count
        upper = sampled_pole_dft_column(p_upper(i), f_hz, sampleRate, sampleCount);
        lower = sampled_pole_dft_column(conj(p_upper(i)), f_hz, sampleRate, sampleCount);
        A(:, i) = upper + lower;
        A(:, pair_count + i) = 1i * (upper - lower);
    end

    for i = 1:real_count
        A(:, 2 * pair_count + i) = sampled_pole_dft_column(p_real(i), f_hz, sampleRate, sampleCount);
    end
end

function G = sampled_pole_dft_column(pole, f_hz, sampleRate, sampleCount)
    a = exp(pole / sampleRate);
    omega = 2 * pi * f_hz(:) / sampleRate;
    q = a * exp(-1i * omega);
    den = 1 - q;
    num = 1 - q .^ sampleCount;
    G = num ./ den;

    nearMask = abs(den) < 1e-10;
    if any(nearMask)
        G(nearMask) = sampleCount;
    end
end

function C = h0_constraint_row(pair_count, real_count)
    C = [2 * ones(1, pair_count), zeros(1, pair_count), ones(1, real_count)];
end

function [r_upper, r_real, p_all, r_all] = unpack_modal_params(theta, p_upper, p_real)
    pair_count = numel(p_upper);
    r_upper = theta(1:pair_count) + 1i * theta(pair_count + (1:pair_count));
    r_real = theta(2 * pair_count + (1:numel(p_real)));

    p_all = [p_upper; conj(p_upper); p_real(:)];
    r_all = [r_upper; conj(r_upper); r_real(:)];
end

function weights = make_relative_weights(H, floorRatio)
    peak = max(abs(H(:)));
    if peak <= eps
        weights = ones(size(H(:)));
        return;
    end

    floorMag = peak * floorRatio;
    weights = 1 ./ max(abs(H(:)), floorMag);
    weights = weights / median(weights);
end

function idx = fit_indices(n, maxFitPoints)
    if maxFitPoints >= n
        idx = (1:n).';
        return;
    end

    idx = unique(round(linspace(1, n, maxFitPoints))).';
    idx(1) = 1;
    idx(end) = n;
end

function idx = plot_indices_for_frequency(f_hz, maxHz, maxPoints)
    all_idx = find(f_hz <= maxHz);
    if isempty(all_idx)
        all_idx = 1;
    end

    if numel(all_idx) > maxPoints
        local = unique(round(linspace(1, numel(all_idx), maxPoints))).';
        idx = all_idx(local);
    else
        idx = all_idx;
    end
end

function y = mag_db(H, cfg)
    y = 20 * log10(max(abs(H(:)), cfg.spectrumMagnitudeFloor));
end

function y = exclusive_cumsum(x)
    x = x(:);
    y = zeros(size(x));
    if numel(x) > 1
        y(2:end) = cumsum(x(1:end-1));
    end
end

function w = one_sided_blackman_window(n, power_value)
    if n <= 1
        w = 0;
        return;
    end

    t = (0:n-1).' / (n - 1);
    u = 0.5 + 0.5 * (t .^ power_value);
    w = blackman_window(u);
    w(end) = 0.0;
end

function w = blackman_window(x)
    w = 0.42 - 0.5 * cos(2.0 * pi * x) + 0.08 * cos(4.0 * pi * x);
    w(x < 0 | x > 1) = 0;
end

function y = apply_dc_compensation(x)
    y = x(:);
    if numel(y) <= 1
        return;
    end

    currentDCSum = sum(y);
    if abs(currentDCSum) < 1e-9
        return;
    end

    window = blackman_harris_window((0:numel(y)-1).' / (numel(y) - 1));
    windowSum = sum(window);
    if abs(windowSum) < 1e-9
        return;
    end

    y = y - window * (currentDCSum / windowSum);
end

function w = blackman_harris_window(x)
    w = 0.35875 ...
        - 0.48829 * cos(2.0 * pi * x) ...
        + 0.14128 * cos(4.0 * pi * x) ...
        - 0.01168 * cos(6.0 * pi * x);
    w(x < 0 | x > 1) = 0;
end

function out = call_without_figures(fun)
    before = findall(groot, 'Type', 'figure');
    oldFigureVisible = get(groot, 'defaultFigureVisible');
    cleanup = onCleanup(@() set(groot, 'defaultFigureVisible', oldFigureVisible));
    set(groot, 'defaultFigureVisible', 'off');

    out = fun();

    after = findall(groot, 'Type', 'figure');
    created = after(~ismember(after, before));
    if ~isempty(created)
        delete(created);
    end

    clear cleanup;
end
