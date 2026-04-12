function result = DesignHybridSystem(n_poles, cfgOverride)
%DESIGNHYBRIDSYSTEM Build FIR tables and fit hybrid IIR residuals.
%
% The default path follows the current HybridBlep design notes:
%   1) Generate min-phase BLIT/BLEP/BLAMP targets.
%   2) Take the first FIR table region and apply the C++ DC compensation.
%   3) Zero-pad FIR tables to the target length.
%   4) Fit FFT(target - zero-padded FIR) with invfreqsTools, then refit residues with
%      fixed poles under h(0+) = 0.

    if nargin < 1 || isempty(n_poles),      n_poles = 8; end
    if nargin < 2 || isempty(cfgOverride),  cfgOverride = struct(); end

    validateattributes(n_poles, {'numeric'}, {'scalar','integer','positive'});

    cfg = default_config();
    cfg = apply_overrides(cfg, cfgOverride);
    cfg = finalize_config(cfg);

    fprintf('====================================================\n');
    fprintf('Hybrid BLIT/BLEP/BLAMP system design\n');
    fprintf('fs / Ts                   = %.12g Hz / %.12g s\n', cfg.fs, cfg.Ts);
    fprintf('n_target                  = %d\n', cfg.n_target);
    fprintf('FIR window/table/total    = %d / %d / %d\n', ...
        cfg.firWindowSize, cfg.firTableNum, cfg.firTotalLength);
    fprintf('bandlimit                 = %.12g\n', cfg.bandlimit);
    fprintf('bandlimit safety factor   = %.12g\n', cfg.bandlimitWc);
    fprintf('bandlimit main-lobe width = %.12g samples\n', cfg.firWindowSize * cfg.bandlimitWc);
    fprintf('target sample rate        = %.12g Hz\n', cfg.targetSampleRate);
    fprintf('IIR pole count            = %d\n', n_poles);
    fprintf('IIR target mode           = %s\n', cfg.iirTargetMode);
    fprintf('Taper IIR target residual = %d\n', logical(cfg.taperIirTargetResidual));
    fprintf('Force IIR target h[0] zero= %d\n', logical(cfg.forceIirTargetInitialZero));
    fprintf('BLAMP dtv                 = %.12g\n', cfg.blampDtv);
    fprintf('====================================================\n\n');

    source = MinphaseBandlimitSiGenerator(cfg.n_target, cfg.bandlimit, cfg.minphaseOversample);
    target = make_target_residuals(source, cfg);

    fir = make_fir_tables(target, cfg);
    fir_plot = plot_fir_table_window(fir.table_saved, cfg);
    prepared = prepare_fit_signals(target, fir, cfg);
    spectra = make_error_spectra(prepared, cfg);
    prefit_plot = plot_prefit_magnitude_spectra(spectra, cfg);
    fits = fit_all_stages(spectra, n_poles, cfg);
    postfit_plot = plot_postfit_magnitude_spectra(spectra, fits, prepared, cfg);

    result = struct();
    result.cfg = cfg;
    result.source = source;
    result.target = target;
    result.fir = fir;
    result.fir_plot = fir_plot;
    result.prepared = prepared;
    result.spectra = spectra;
    result.prefit_plot = prefit_plot;
    result.fits = fits;
    result.postfit_plot = postfit_plot;
end

function cfg = default_config()
    cfg = struct();
    cfg.fs = 48000;
    cfg.Ts = [];
    cfg.n_target = 131072 * 8;
    cfg.firWindowSize = 16;
    cfg.firTableNum = 64;
    cfg.firTotalLength = [];
    cfg.bandlimit = [];
    cfg.bandlimitWc = 0.95;
    cfg.minphaseOversample = 8;
    cfg.causalWindowPower = 10.0;
    cfg.targetSampleRate = [];
    cfg.blampDtv = [];
    cfg.maxFitPoints = 16384;
    cfg.relativeFitWeightFloor = 1e-4;
    cfg.hideInvfreqsToolFigures = true;
    cfg.runInvfreqsFit = true;
    cfg.spectrumPlotMaxHz = [];
    cfg.spectrumPlotPoints = 4096;
    cfg.spectrumMagnitudeFloor = 1e-300;
    cfg.iirImpulsePlotSamples = [];
    cfg.enforceStablePoles = true;
    cfg.minPoleDampingRad = [];
    cfg.iirTargetMode = 'parallel-residual';
    cfg.taperIirTargetResidual = false;
    cfg.forceIirTargetInitialZero = false;
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
    validateattributes(cfg.n_target, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.firWindowSize, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.firTableNum, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.minphaseOversample, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(cfg.causalWindowPower, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.relativeFitWeightFloor, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.bandlimitWc, {'numeric'}, {'scalar','positive'});
    validateattributes(cfg.spectrumPlotPoints, {'numeric'}, {'scalar','integer','>=',2});
    validateattributes(cfg.spectrumMagnitudeFloor, {'numeric'}, {'scalar','positive'});
    if mod(cfg.n_target, 2) ~= 0
        error('n_target must be even so the one-sided FFT spectrum includes an exact Nyquist bin.');
    end

    cfg.Ts = 1 / cfg.fs;
    cfg.firTotalLength = cfg.firWindowSize * cfg.firTableNum;
    if cfg.firTotalLength > cfg.n_target
        error('firTotalLength (%d) must not exceed n_target (%d).', cfg.firTotalLength, cfg.n_target);
    end

    if isempty(cfg.bandlimit)
        cfg.bandlimit = 2 / (cfg.firWindowSize * cfg.bandlimitWc);
    end
    validateattributes(cfg.bandlimit, {'numeric'}, {'scalar','positive','<',1});

    if isempty(cfg.targetSampleRate)
        cfg.targetSampleRate = cfg.fs * cfg.firTableNum;
    end
    validateattributes(cfg.targetSampleRate, {'numeric'}, {'scalar','positive'});
    if isempty(cfg.spectrumPlotMaxHz)
        cfg.spectrumPlotMaxHz = min(cfg.targetSampleRate / 2, cfg.bandlimit * cfg.targetSampleRate);
    end
    validateattributes(cfg.spectrumPlotMaxHz, {'numeric'}, {'scalar','positive'});

    if isempty(cfg.blampDtv)
        cfg.blampDtv = 1 / cfg.firTableNum;
    end
    validateattributes(cfg.blampDtv, {'numeric'}, {'scalar','positive'});
    if isempty(cfg.iirImpulsePlotSamples)
        cfg.iirImpulsePlotSamples = min(cfg.n_target, max(cfg.firTotalLength, 2048));
    end
    validateattributes(cfg.iirImpulsePlotSamples, {'numeric'}, {'scalar','integer','positive'});
    if isempty(cfg.minPoleDampingRad)
        cfg.minPoleDampingRad = cfg.targetSampleRate * 1e-6;
    end
    validateattributes(cfg.minPoleDampingRad, {'numeric'}, {'scalar','positive'});

    if ~(islogical(cfg.hideInvfreqsToolFigures) || ...
            (isnumeric(cfg.hideInvfreqsToolFigures) && isscalar(cfg.hideInvfreqsToolFigures)))
        error('hideInvfreqsToolFigures must be a scalar logical/numeric flag.');
    end
    if ~(islogical(cfg.runInvfreqsFit) || ...
            (isnumeric(cfg.runInvfreqsFit) && isscalar(cfg.runInvfreqsFit)))
        error('runInvfreqsFit must be a scalar logical/numeric flag.');
    end
    if ~(islogical(cfg.enforceStablePoles) || ...
            (isnumeric(cfg.enforceStablePoles) && isscalar(cfg.enforceStablePoles)))
        error('enforceStablePoles must be a scalar logical/numeric flag.');
    end
    if ~(ischar(cfg.iirTargetMode) || isstring(cfg.iirTargetMode))
        error('iirTargetMode must be a character vector or string scalar.');
    end
    cfg.iirTargetMode = char(cfg.iirTargetMode);
    if ~(islogical(cfg.taperIirTargetResidual) || ...
            (isnumeric(cfg.taperIirTargetResidual) && isscalar(cfg.taperIirTargetResidual)))
        error('taperIirTargetResidual must be a scalar logical/numeric flag.');
    end
    if ~(islogical(cfg.forceIirTargetInitialZero) || ...
            (isnumeric(cfg.forceIirTargetInitialZero) && isscalar(cfg.forceIirTargetInitialZero)))
        error('forceIirTargetInitialZero must be a scalar logical/numeric flag.');
    end
end

function target = make_target_residuals(source, cfg)
    sinc_unit_area = source.minphase_sinc(:);

    target = struct();
    target.names = {'blit', 'blep', 'blamp'};
    target.blit = sinc_unit_area * cfg.firTableNum;
    target.blep = cumsum(sinc_unit_area) - 1.0;

    % Integrate from the sample boundary so BLAMP starts at 0.
    target.blamp = exclusive_cumsum(target.blep) * cfg.blampDtv;
end

function fir = make_fir_tables(target, cfg)
    fir = struct();
    fir.raw = struct();
    fir.dc_compensated = struct();

    for i = 1:numel(target.names)
        name = target.names{i};
        raw = target.(name)(1:cfg.firTotalLength);
        fir.raw.(name) = raw;
        fir.dc_compensated.(name) = apply_dc_compensation(raw);
    end

    fir.table_saved = fir.dc_compensated;
end

function plot_data = plot_fir_table_window(firTable, cfg)
    x = (0:cfg.firTotalLength-1).' / cfg.firTableNum;
    names = {'blit', 'blep', 'blamp'};

    figure('Name', 'Hybrid FIR table over full firWindowSize');
    for i = 1:numel(names)
        name = names{i};
        subplot(numel(names), 1, i);
        plot(x, firTable.(name), 'LineWidth', 1.2);
        grid on;
        xlabel('Sample offset');
        ylabel(name);
        title(sprintf('FIR %s table, full %.12g-sample window', upper(name), cfg.firWindowSize));
        xlim([0, cfg.firWindowSize]);
    end

    plot_data = struct();
    plot_data.x_sample_offset = x;
    plot_data.blit = firTable.blit;
    plot_data.blep = firTable.blep;
    plot_data.blamp = firTable.blamp;
end

function prepared = prepare_fit_signals(target, fir, cfg)
    prepared = struct();
    prepared.target = struct();
    prepared.fir = struct();
    prepared.error_time = struct();
    prepared.error_time_raw = struct();

    for i = 1:numel(target.names)
        name = target.names{i};
        switch lower(cfg.iirTargetMode)
            case {'parallel-residual', 'tail'}
                target_prepared = target.(name)(:);
                fir_prepared = pad_with_value(fir.dc_compensated.(name), cfg.n_target, 0.0);
                error_raw = target_prepared - fir_prepared;

                error_time = error_raw;
                if cfg.forceIirTargetInitialZero
                    error_time(1) = 0.0;
                end
                if cfg.taperIirTargetResidual
                    error_time = error_time .* one_sided_blackman_window(numel(error_time), cfg.causalWindowPower);
                end

            case 'legacy-zero-start'
                target_prepared = zero_start_and_taper(target.(name), cfg);
                fir_padded = pad_with_value(fir.dc_compensated.(name), cfg.n_target, fir.dc_compensated.(name)(end));
                fir_prepared = zero_start_and_taper(fir_padded, cfg);
                error_raw = target_prepared - fir_prepared;
                error_time = error_raw;

            otherwise
                error('Unknown iirTargetMode "%s".', cfg.iirTargetMode);
        end

        prepared.target.(name) = target_prepared;
        prepared.fir.(name) = fir_prepared;
        prepared.error_time_raw.(name) = error_raw;
        prepared.error_time.(name) = error_time;
    end
end

function spectra = make_error_spectra(prepared, cfg)
    names = fieldnames(prepared.error_time);
    spectra = struct();
    spectra.freq_hz_full = (0:cfg.n_target/2).' * (cfg.targetSampleRate / cfg.n_target);
    spectra.freq_hz_fit = struct();
    spectra.target_fft = struct();
    spectra.fir_fft = struct();
    spectra.error_fft = struct();
    spectra.error_fit = struct();
    spectra.fit_indices = struct();

    for i = 1:numel(names)
        name = names{i};
        target_fft = fft(prepared.target.(name));
        fir_fft = fft(prepared.fir.(name));
        error_fft = target_fft - fir_fft;

        half = cfg.n_target / 2 + 1;
        spectra.target_fft.(name) = target_fft(1:half);
        spectra.fir_fft.(name) = fir_fft(1:half);
        spectra.error_fft.(name) = error_fft(1:half);

        idx = fit_indices(half, cfg.maxFitPoints);
        spectra.fit_indices.(name) = idx;
        spectra.freq_hz_fit.(name) = spectra.freq_hz_full(idx);
        spectra.error_fit.(name) = spectra.error_fft.(name)(idx);
    end
end

function plot_data = plot_prefit_magnitude_spectra(spectra, cfg)
    names = fieldnames(spectra.target_fft);
    f_hz = spectra.freq_hz_full;
    freqMask = f_hz <= cfg.spectrumPlotMaxHz;

    figure('Name', 'Hybrid prefit FIR vs target magnitude spectra');
    for i = 1:numel(names)
        name = names{i};
        target_fft = spectra.target_fft.(name);
        fir_fft = spectra.fir_fft.(name);
        lastIdx = find(freqMask, 1, 'last');

        subplot(numel(names), 1, i);
        plot(f_hz(freqMask), mag_db(target_fft(freqMask), cfg), 'LineWidth', 1.2);
        hold on;
        plot(f_hz(freqMask), mag_db(fir_fft(freqMask), cfg), '--', 'LineWidth', 1.1);
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('%s target vs short-window FIR before IIR fit', upper(name)));
        legend({'target', 'short FIR'}, 'Location', 'best');
        xlim([f_hz(1), f_hz(lastIdx)]);
    end

    plot_data = struct();
    plot_data.freq_hz = f_hz;
    plot_data.freq_mask = freqMask;
end

function plot_data = plot_postfit_magnitude_spectra(spectra, fits, prepared, cfg)
    names = fieldnames(spectra.target_fft);
    f_hz = spectra.freq_hz_full;
    freqMask = f_hz <= cfg.spectrumPlotMaxHz;

    plot_data = struct();
    plot_data.freq_hz = f_hz;
    plot_data.freq_mask = freqMask;
    plot_data.iir_fft_direct = struct();
    plot_data.iir_fft_sampled = struct();
    plot_data.iir_time = struct();
    plot_data.hybrid_time = struct();
    plot_data.hybrid_fft_direct = struct();
    plot_data.hybrid_fft_sampled = struct();
    plot_data.iir_target_time = struct();
    plot_data.valid = false;

    if ~has_any_valid_h0_fit(fits, names)
        fprintf('Skipping postfit FIR+IIR magnitude plot because no valid IIR fit is available.\n');
        return;
    end

    impulseFigure = figure('Name', 'Hybrid IIR residual target vs sampled IIR');
    postfitFigure = figure('Name', 'Hybrid postfit FIR+IIR vs target magnitude spectra');
    for i = 1:numel(names)
        name = names{i};
        target_fft = spectra.target_fft.(name);
        fir_fft = spectra.fir_fft.(name);
        lastIdx = find(freqMask, 1, 'last');

        figure(postfitFigure);
        subplot(numel(names), 1, i);
        plot(f_hz(freqMask), mag_db(target_fft(freqMask), cfg), 'LineWidth', 1.2);
        hold on;

        if isfield(fits, name) && isfield(fits.(name), 'h0_zero') && ...
                isfield(fits.(name).h0_zero, 'valid') && fits.(name).h0_zero.valid
            iir_time = sampled_modal_time_response(fits.(name).h0_zero, cfg.n_target, cfg.targetSampleRate);
            iir_fft_sampled_full = fft(iir_time);
            iir_fft_sampled = iir_fft_sampled_full(1:numel(f_hz));
            hybrid_time = prepared.fir.(name) + iir_time;
            hybrid_fft_full = fft(hybrid_time);
            hybrid_fft_sampled = hybrid_fft_full(1:numel(f_hz));
            iir_fft_direct = sampled_modal_frequency_response(fits.(name).h0_zero, f_hz, cfg);
            hybrid_fft_direct = fir_fft + iir_fft_direct;

            plot(f_hz(freqMask), mag_db(fir_fft(freqMask), cfg), ':', 'LineWidth', 1.0);
            plot(f_hz(freqMask), mag_db(hybrid_fft_sampled(freqMask), cfg), ...
                '--', 'LineWidth', 1.1);
            plot_data.iir_fft_direct.(name) = iir_fft_direct;
            plot_data.iir_fft_sampled.(name) = iir_fft_sampled;
            plot_data.iir_time.(name) = iir_time;
            plot_data.hybrid_time.(name) = hybrid_time;
            plot_data.hybrid_fft_direct.(name) = hybrid_fft_direct;
            plot_data.hybrid_fft_sampled.(name) = hybrid_fft_sampled;
            plot_data.iir_target_time.(name) = prepared.error_time.(name);
            legendEntries = {'target', 'short FIR', 'short FIR + sampled IIR'};

            figure(impulseFigure);
            subplot(numel(names), 1, i);
            iirPlotSamples = min(cfg.iirImpulsePlotSamples, numel(iir_time));
            x = (0:iirPlotSamples-1).' / cfg.firTableNum;
            plot(x, prepared.error_time.(name)(1:iirPlotSamples), 'LineWidth', 1.1);
            hold on;
            plot(x, iir_time(1:iirPlotSamples), '--', 'LineWidth', 1.1);
            grid on;
            xlabel('Sample offset');
            ylabel(name);
            title(sprintf('%s IIR residual target vs sampled IIR', upper(name)));
            legend({'IIR target', 'sampled IIR'}, 'Location', 'best');
        else
            plot(f_hz(freqMask), mag_db(fir_fft(freqMask), cfg), '--', 'LineWidth', 1.1);
            plot_data.iir_fft_direct.(name) = [];
            plot_data.iir_fft_sampled.(name) = [];
            plot_data.iir_time.(name) = [];
            plot_data.hybrid_time.(name) = [];
            plot_data.hybrid_fft_direct.(name) = [];
            plot_data.hybrid_fft_sampled.(name) = [];
            plot_data.iir_target_time.(name) = [];
            legendEntries = {'target', 'short FIR only'};
        end

        figure(postfitFigure);
        subplot(numel(names), 1, i);
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('%s target vs hybrid FIR+sampled IIR after IIR fit', upper(name)));
        legend(legendEntries, 'Location', 'best');
        xlim([f_hz(1), f_hz(lastIdx)]);
    end

    plot_data.valid = true;
end

function tf = has_any_valid_h0_fit(fits, names)
    tf = false;
    for i = 1:numel(names)
        name = names{i};
        if isfield(fits, name) && isfield(fits.(name), 'h0_zero') && ...
                isfield(fits.(name).h0_zero, 'valid') && fits.(name).h0_zero.valid
            tf = true;
            return;
        end
    end
end

function H = sampled_modal_frequency_response(model, f_hz, cfg)
    A = sampled_modal_basis_matrix(f_hz, model.upper_poles_rad(:), ...
        model.real_poles_rad(:), cfg.targetSampleRate, cfg.n_target);
    theta = pack_modal_params(model.upper_residues(:), model.real_residues(:));
    H = A * theta;
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

function y = mag_db(H, cfg)
    y = 20 * log10(max(abs(H(:)), cfg.spectrumMagnitudeFloor));
end

function idx = fit_indices(n, maxFitPoints)
    if ~isfinite(maxFitPoints) || maxFitPoints >= n
        idx = (1:n).';
        return;
    end

    validateattributes(maxFitPoints, {'numeric'}, {'scalar','integer','>=',2});
    idx = unique(round(linspace(1, n, maxFitPoints))).';
    idx(1) = 1;
    idx(end) = n;
end

function fits = fit_all_stages(spectra, n_poles, cfg)
    names = fieldnames(spectra.error_fit);
    fits = struct();

    if ~cfg.runInvfreqsFit
        fprintf('Skipping invfreqsTools fitting because cfg.runInvfreqsFit is false.\n');
        for i = 1:numel(names)
            name = names{i};
            fits.(name) = struct();
            fits.(name).tool = struct();
            fits.(name).h0_zero = struct();
        end
        return;
    end

    for i = 1:numel(names)
        name = names{i};
        fprintf('Fitting %s error spectrum with invfreqsTools...\n', upper(name));
        H = spectra.error_fit.(name);
        f_hz = spectra.freq_hz_fit.(name);
        tool_fit = call_invfreqs_tools(H, n_poles, cfg.targetSampleRate, f_hz, cfg);

        constrained = struct();
        if isfield(tool_fit, 'invfreqs') && tool_fit.invfreqs.valid
            constrained = refit_residues_h0_zero(tool_fit.invfreqs, H, f_hz, cfg);
        else
            warning('%s invfreqsTools fit did not return a valid invfreqs model.', upper(name));
        end

        fits.(name) = struct();
        fits.(name).tool = tool_fit;
        fits.(name).h0_zero = constrained;
    end
end

function fit = call_invfreqs_tools(H, n_poles, sample_rate, f_hz, cfg)
    oldFigureVisible = get(groot, 'defaultFigureVisible');
    cleanup = onCleanup(@() set(groot, 'defaultFigureVisible', oldFigureVisible));
    if cfg.hideInvfreqsToolFigures
        set(groot, 'defaultFigureVisible', 'off');
    end

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

function constrained = refit_residues_h0_zero(model, H, f_hz, cfg)
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
    Z = null(C);
    if isempty(Z)
        theta = zeros(size(C, 2), 1);
    else
        theta = Z * ((Ar * Z) \ br);
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

function theta = pack_modal_params(r_upper, r_real)
    r_upper = r_upper(:);
    r_real = r_real(:);
    theta = [real(r_upper); imag(r_upper); real(r_real)];
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

function y = pad_with_value(x, n, value)
    x = x(:);
    if numel(x) > n
        error('Cannot pad signal of length %d down to %d.', numel(x), n);
    end

    y = zeros(n, 1);
    y(1:numel(x)) = x;
    if numel(x) < n
        y(numel(x) + 1:end) = value;
    end
end

function y = zero_start_and_taper(x, cfg)
    y = x(:) - x(1);
    y = y .* one_sided_blackman_window(numel(y), cfg.causalWindowPower);
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
