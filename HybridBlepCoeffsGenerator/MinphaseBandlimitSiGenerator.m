function result = MinphaseBandlimitSiGenerator(n, bandlimit, oversample)
%MINPHASEBANDLIMITSIGENERATOR Generate bandlimited min-phase sinc and SI.
%
% Default settings:
%   n         = 32768
%   bandlimit = 0.495
%   oversample = 8

    if nargin < 1 || isempty(n),          n = 32768; end
    if nargin < 2 || isempty(bandlimit),  bandlimit = 0.495; end
    if nargin < 3 || isempty(oversample), oversample = 8; end

    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(bandlimit, {'numeric'}, {'scalar','positive','<',1});
    validateattributes(oversample, {'numeric'}, {'scalar','integer','positive'});
    if mod(n * oversample, 2) ~= 0
        error('n * oversample must be even for the real-cepstrum causal window.');
    end

    cfg = struct();
    cfg.n = n;
    cfg.bandlimit = bandlimit;
    cfg.oversample = oversample;
    cfg.fft_len = n * oversample;
    cfg.causal_window_power = 10.0;
    cfg.plot_samples = min(100, n);

    i = (0:n-1).';
    x = i - floor(n / 2);
    xw = i / (n - 1);

    seed_window = blackman_window(xw);
    seed_sinc = sinc_arg(pi * bandlimit * x) .* seed_window;

    minphase_sinc_full = minimum_phase_from_padded_seed(seed_sinc, cfg.fft_len);
    minphase_sinc_raw = minphase_sinc_full(1:n);

    t = i / (n - 1);
    causal_window = one_sided_blackman_window(t, cfg.causal_window_power);
    minphase_sinc_windowed = minphase_sinc_raw .* causal_window;

    sinc_area = sum(minphase_sinc_windowed);
    if abs(sinc_area) < 1e-18
        error('Windowed min-phase sinc has too little DC area to normalize.');
    end
    minphase_sinc = minphase_sinc_windowed / sinc_area;

    minphase_si_step = cumsum(minphase_sinc);
    minphase_si = minphase_si_step - 1.0;
    minphase_int_si = exclusive_cumsum(minphase_si);

    fprintf('====================================================\n');
    fprintf('Min-phase bandlimited sinc/SI generator\n');
    fprintf('n                         = %d\n', n);
    fprintf('bandlimit                 = %.12g\n', bandlimit);
    fprintf('x range                   = %.12g .. %.12g\n', x(1), x(end));
    fprintf('cepstrum oversample       = %d x\n', oversample);
    fprintf('cepstrum fft length       = %d\n', cfg.fft_len);
    fprintf('causal window power       = %.12g\n', cfg.causal_window_power);
    fprintf('raw sinc sum              = %.12g\n', sum(seed_sinc));
    fprintf('windowed minphase sum     = %.12g\n', sinc_area);
    fprintf('normalized sinc sum       = %.12g\n', sum(minphase_sinc));
    fprintf('final SI step value       = %.12g\n', minphase_si_step(end));
    fprintf('final SI residual value   = %.12g\n', minphase_si(end));
    fprintf('final int SI value        = %.12g\n', minphase_int_si(end));
    fprintf('causal window first/last  = %.12g / %.12g\n', causal_window(1), causal_window(end));
    fprintf('====================================================\n\n');

    plot_first_samples(minphase_sinc, minphase_si, minphase_int_si, cfg);

    result = struct();
    result.cfg = cfg;
    result.x = x;
    result.seed_window = seed_window;
    result.seed_sinc = seed_sinc;
    result.minphase_sinc_full = minphase_sinc_full;
    result.minphase_sinc_raw = minphase_sinc_raw;
    result.causal_window = causal_window;
    result.minphase_sinc_windowed = minphase_sinc_windowed;
    result.minphase_sinc = minphase_sinc;
    result.minphase_si_step = minphase_si_step;
    result.minphase_si = minphase_si;
    result.minphase_int_si = minphase_int_si;
end

function y = sinc_arg(x)
    y = ones(size(x));
    mask = abs(x) >= 1e-12;
    y(mask) = sin(x(mask)) ./ x(mask);
end

function y = minimum_phase_from_padded_seed(seed, fft_len)
    x = zeros(fft_len, 1);
    x(1:numel(seed)) = seed(:);

    X = fft(x);
    cep = ifft(log(abs(X) + 1e-100), 'symmetric');

    % Causal real-cepstrum lifter:
    % keep cep[0], double positive quefrencies, keep Nyquist, zero the rest.
    cep(2:fft_len/2) = cep(2:fft_len/2) * 2.0;
    cep(fft_len/2 + 2:end) = 0.0;

    Y = exp(fft(cep));
    y = real(ifft(Y));
end

function w = one_sided_blackman_window(t, power_value)
    u = 0.5 + 0.5 * (t .^ power_value);
    w = blackman_window(u);
    w(end) = 0.0;
end

function w = blackman_window(x)
    w = 0.42 - 0.5 * cos(2.0 * pi * x) + 0.08 * cos(4.0 * pi * x);
    w(x < 0 | x > 1) = 0;
end

function y = exclusive_cumsum(x)
    y = zeros(size(x));
    if numel(x) > 1
        y(2:end) = cumsum(x(1:end-1));
    end
end

function plot_first_samples(minphase_sinc, minphase_si, minphase_int_si, cfg)
    idx = (1:cfg.plot_samples).';
    sample = idx - 1;

    figure('Name', 'Minphase bandlimited sinc/SI/int SI first 100 samples');

    subplot(3, 1, 1);
    plot(sample, minphase_sinc(idx), 'LineWidth', 1.2);
    grid on;
    xlabel('Sample');
    ylabel('h');
    title('Min-phase bandlimited sinc');

    subplot(3, 1, 2);
    plot(sample, minphase_si(idx), 'LineWidth', 1.2);
    grid on;
    xlabel('Sample');
    ylabel('SI');
    title('Min-phase SI residual');

    subplot(3, 1, 3);
    plot(sample, minphase_int_si(idx), 'LineWidth', 1.2);
    grid on;
    xlabel('Sample');
    ylabel('int SI');
    title('Integrated min-phase SI');
end
