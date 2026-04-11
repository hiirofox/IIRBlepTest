function result = blampresiduecoeffsgenerator(n, fc_hz, Rp, Rs, fs, dc_block_rad, dc_block_order, dc_spread_oct, response_samples)
% Analog BLAMP-residual modal export for direct InjectImpulse() playback.
%
% Transfer function:
%   R(s) = fs * (Hlp(s) - 1) / s^2 * Hantidc(s)
%
% The extra factor fs converts the analog time integral into a sample-domain
% BLAMP residual. Without it, the response is smaller by Ts and looks like 0
% in the C++ plots.
%
% Defaults:
%   Hlp      : 12th-order analog elliptic low-pass, DC-normalized to Hlp(0)=1
%   Hhp      : 12th-order analog elliptic high-pass, generated for comparison
%   Hantidc : 2nd-order DC blocker, prod_k s / (s + wd_k)
%
% Export format:
%   twoPoleParams = [pre, pim(rad/s), rre, rim, ...]
%   onePoleParams = [pre, rre, ...]
%
% Intended C++ use:
%   feed the exported residues to InjectImpulse() to generate the BLAMP residual.
%   If residue() returns a direct term k_dir for a future formula variant, apply
%   that direct term at the injection sample in addition to the modal tail.

    if nargin < 1 || isempty(n),                 n = 12; end
    if nargin < 2 || isempty(fc_hz),             fc_hz = 23500; end
    if nargin < 3 || isempty(Rp),                Rp = 1; end
    if nargin < 4 || isempty(Rs),                Rs = 120; end
    if nargin < 5 || isempty(fs),                fs = 48000; end
    if nargin < 6 || isempty(dc_block_rad),      dc_block_rad = 400.0; end
    if nargin < 7 || isempty(dc_block_order),    dc_block_order = 2; end
    if nargin < 8 || isempty(dc_spread_oct),     dc_spread_oct = 1; end
    if nargin < 9 || isempty(response_samples),  response_samples = 500; end

    validateattributes(n,                {'numeric'}, {'scalar','integer','positive'});
    validateattributes(fc_hz,            {'numeric'}, {'scalar','positive'});
    validateattributes(Rp,               {'numeric'}, {'scalar','positive'});
    validateattributes(Rs,               {'numeric'}, {'scalar','positive'});
    validateattributes(fs,               {'numeric'}, {'scalar','positive'});
    validateattributes(dc_block_order,   {'numeric'}, {'scalar','integer','positive'});
    validateattributes(dc_spread_oct,    {'numeric'}, {'scalar','nonnegative'});
    validateattributes(response_samples, {'numeric'}, {'scalar','integer','positive'});

    if isscalar(dc_block_rad)
        validateattributes(dc_block_rad, {'numeric'}, {'scalar','positive'});
    else
        validateattributes(dc_block_rad, {'numeric'}, {'vector','positive'});
    end

    wc = 2*pi*fc_hz;

    [z, p0, k0] = ellipap(n, Rp, Rs);
    [b0, a0] = zp2tf(z, p0, k0);
    [b_lp, a_lp] = lp2lp(b0, a0, wc);
    [b_hp, a_hp] = lp2hp(b0, a0, wc);

    Hlp0_before = polyval(b_lp, 0) / polyval(a_lp, 0);
    b_lp = b_lp / Hlp0_before;
    Hlp0 = polyval(b_lp, 0) / polyval(a_lp, 0);

    HhpInf_before = high_frequency_gain(b_hp, a_hp);
    if ~isfinite(HhpInf_before) || abs(HhpInf_before) <= eps
        error('Cannot normalize Hhp(inf); high-frequency gain is %.12g %+.12gj.', real(HhpInf_before), imag(HhpInf_before));
    end
    b_hp = b_hp / HhpInf_before;
    HhpInf = high_frequency_gain(b_hp, a_hp);
    Hhp0 = polyval(b_hp, 0) / polyval(a_hp, 0);

    dc_poles_rad = build_dc_poles(dc_block_rad, dc_block_order, dc_spread_oct);
    [b_dc, a_dc] = make_dc_blocker(dc_poles_rad);

    % Direct transfer-function expression. This is intentionally written in
    % the same form as the design idea so the sign and the direct term are clear.
    s = tf('s');
    Hlp = tf(b_lp, a_lp);
    Hhp = tf(b_hp, a_hp);
    Hantidc = tf(b_dc, a_dc);
    % Alternative explicit high-pass residual:
    % R = minreal((Hlp-1)/s*Hantidc, 1e-7);
    R = minreal(fs*(Hlp-1)/(s*s)*Hantidc, 1e-7);
    [b, a] = tfdata(R, 'v');
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);

    area = polyval(b, 0) / polyval(a, 0);

    [r_all, p_all, k_dir] = residue(b, a);
    directGain = sum(k_dir);
    h0_modal = real(sum(r_all));
    h0_display = directGain + h0_modal;
    slope0_per_second = initial_derivative(b, a);
    slope0_per_sample = slope0_per_second / fs;
    slope0_error_from_minus_one = slope0_per_sample + 1;
    exportResidueNorm = abs(slope0_per_sample);
    if ~isfinite(exportResidueNorm) || exportResidueNorm <= eps
        error('Cannot normalize exported residues; sample-domain initial slope is %.12g %+.12gj.', real(slope0_per_sample), imag(slope0_per_sample));
    end
    directGainExport = directGain / exportResidueNorm;
    h0_export = h0_display / exportResidueNorm;
    slope0_export = slope0_per_sample / exportResidueNorm;
    slope0_export_error_from_minus_one = slope0_export + 1;
    [p_two, r_two, p_one, r_one] = split_modal_poles(p_all, r_all);
    r_two_export = r_two / exportResidueNorm;
    r_one_export = r_one / exportResidueNorm;

    twoPoleParams = zeros(numel(p_two) * 4, 1);
    for i = 1:numel(p_two)
        twoPoleParams((i-1)*4 + 1) = real(p_two(i));
        twoPoleParams((i-1)*4 + 2) = imag(p_two(i));
        twoPoleParams((i-1)*4 + 3) = real(r_two_export(i));
        twoPoleParams((i-1)*4 + 4) = imag(r_two_export(i));
    end

    onePoleParams = zeros(numel(p_one) * 2, 1);
    for i = 1:numel(p_one)
        onePoleParams((i-1)*2 + 1) = p_one(i);
        onePoleParams((i-1)*2 + 2) = real(r_one_export(i));
    end

    t_samples = (0:response_samples-1) / fs;
    r_all_export = r_all / exportResidueNorm;
    h_samples = real(sum(bsxfun(@times, r_all_export(:), exp(p_all(:) * t_samples)), 1));
    h_samples(1) = h_samples(1) + directGainExport;
    t_analog = linspace(0, (response_samples-1) / fs, max(4096, response_samples * 8));
    h_analog = real(sum(bsxfun(@times, r_all_export(:), exp(p_all(:) * t_analog)), 1));
    h_analog(1) = h_analog(1) + directGainExport;

    f_hz = unique([0, logspace(-3, log10(max(10 * fs / 2, fc_hz * 2)), 16384)]);
    H = freqs(b, a, 2*pi*f_hz);
    H = H(:).';

    % BLAMP residual is already zero-initial. Keep DC out of the log plot and
    % show the transfer function directly.
    H_plot = NaN(size(H));
    positive_freq_mask = f_hz > 0;
    H_plot(positive_freq_mask) = H(positive_freq_mask);
    H_abs = abs(H_plot);
    H_peak = max(H_abs(positive_freq_mask));
    H_db_rel = 20*log10(max(H_abs / H_peak, 1e-300));
    H_nyquist = freqs(b, a, [2*pi*(fs/2), 2*pi*(fs/2)]);
    H_nyquist = H_nyquist(1);
    H_nyquist_plot = H_nyquist;
    H_nyquist_db_rel = 20*log10(max(abs(H_nyquist_plot) / H_peak, 1e-300));

    fprintf('====================================================\n');
    fprintf('BLAMP residual modal export\n');
    fprintf('Formula                   = fs * (Hlp(s) - 1) / s^2 * Hantidc(s)\n');
    fprintf('Elliptic order n          = %d\n', n);
    fprintf('Cutoff fc                 = %.12g Hz\n', fc_hz);
    fprintf('Passband ripple Rp        = %.12g dB\n', Rp);
    fprintf('Stopband attenuation Rs   = %.12g dB\n', Rs);
    fprintf('Hlp(0) before normalize   = %.12g %+.12gj\n', real(Hlp0_before), imag(Hlp0_before));
    fprintf('Hlp(0) after normalize    = %.12g %+.12gj\n', real(Hlp0), imag(Hlp0));
    fprintf('Hhp(inf) before normalize = %.12g %+.12gj\n', real(HhpInf_before), imag(HhpInf_before));
    fprintf('Hhp(inf) after normalize  = %.12g %+.12gj\n', real(HhpInf), imag(HhpInf));
    fprintf('Hhp(0)                    = %.12g %+.12gj\n', real(Hhp0), imag(Hhp0));
    fprintf('DC blocker poles (rad/s)  = ');
    fprintf('%.12g ', dc_poles_rad);
    fprintf('\n');
    fprintf('Frequency plot max        = %.12g Hz (10 * Nyquist)\n', 10 * fs / 2);
    fprintf('Frequency plot view       = R(s)\n');
    fprintf('Nyquist suppression       = %.12g dB relative to response peak\n', H_nyquist_db_rel);
    fprintf('directGain k_dir          = %.12g %+.12gj\n', real(directGain), imag(directGain));
    fprintf('modal tail h(0+)          = %.12g\n', h0_modal);
    fprintf('display h[0] direct+tail  = %.12g %+.12gj\n', real(h0_display), imag(h0_display));
    fprintf('initial slope/sample      = %.12g %+.12gj\n', real(slope0_per_sample), imag(slope0_per_sample));
    fprintf('slope/sample - (-1)       = %.12g %+.12gj\n', real(slope0_error_from_minus_one), imag(slope0_error_from_minus_one));
    fprintf('C++ residue export norm   = abs(slope/sample) = %.12g\n', exportResidueNorm);
    fprintf('C++ export h(0+)          = %.12g %+.12gj\n', real(h0_export), imag(h0_export));
    fprintf('C++ export slope/sample   = %.12g %+.12gj\n', real(slope0_export), imag(slope0_export));
    fprintf('C++ export slope - (-1)   = %.12g %+.12gj\n', real(slope0_export_error_from_minus_one), imag(slope0_export_error_from_minus_one));
    fprintf('integral h(t) dt = H(0)  = %.12g %+.12gj\n', real(area), imag(area));
    fprintf('Two-pole modal count      = %d\n', numel(p_two));
    fprintf('One-pole modal count      = %d\n', numel(p_one));
    fprintf('====================================================\n\n');

    print_cpp_vector(twoPoleParams, 'twoPoleParams', 4);
    fprintf('\n');
    print_cpp_vector(onePoleParams, 'onePoleParams', 2);
    fprintf('\nfloat directGain = %.12gf;\n', real(directGainExport));

    if any(abs(k_dir) > 1e-9)
        fprintf('\n%% NOTE: residue() returned direct term k_dir:\n');
        disp(k_dir);
    end

    figure('Name', 'BLAMP residual modal export');

    subplot(2, 1, 1);
    semilogx(f_hz(2:end), H_db_rel(2:end), 'LineWidth', 1.2);
    hold on;
    yl = ylim;
    plot([fs/2, fs/2], yl, '--', 'LineWidth', 1.1);
    ylim(yl);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB, peak-normalized)');
    title(sprintf('R(s), Nyquist = %.2f dB', H_nyquist_db_rel));
    xlim([f_hz(2), f_hz(end)]);

    subplot(2, 1, 2);
    plot(t_analog * 1000, h_analog, 'LineWidth', 1.2);
    hold on;
    plot([t_analog(1), t_analog(end)] * 1000, [0, 0], ':');
    plot([t_analog(1), t_analog(end)] * 1000, [-1, -1], '--');
    plot(t_samples * 1000, h_samples, '.', 'MarkerSize', 6);
    grid on;
    xlabel('Time (ms)');
    ylabel('h(t)');
    title(sprintf('Sample-domain BLAMP residual, time span = %d samples @ %.0f Hz', response_samples, fs));
    xlim([0, t_analog(end) * 1000]);

    result = struct();
    result.b = b;
    result.a = a;
    result.b_lp = b_lp;
    result.a_lp = a_lp;
    result.b_hp = b_hp;
    result.a_hp = a_hp;
    result.b_dc = b_dc;
    result.a_dc = a_dc;
    result.poles_all_rad = p_all;
    result.residues_all = r_all;
    result.k_dir = k_dir;
    result.two_poles_rad = p_two;
    result.two_residues = r_two;
    result.one_poles_rad = p_one;
    result.one_residues = r_one;
    result.two_residues_export = r_two_export;
    result.one_residues_export = r_one_export;
    result.twoPoleParams = twoPoleParams;
    result.onePoleParams = onePoleParams;
    result.dc_block_poles_rad = dc_poles_rad;
    result.f_hz = f_hz;
    result.H = H;
    result.H_plot = H_plot;
    result.H_db_rel = H_db_rel;
    result.H_peak = H_peak;
    result.H_nyquist = H_nyquist;
    result.H_nyquist_plot = H_nyquist_plot;
    result.H_nyquist_db_rel = H_nyquist_db_rel;
    result.t_samples = t_samples;
    result.h_samples = h_samples;
    result.t_analog = t_analog;
    result.h_analog = h_analog;
    result.directGain = directGain;
    result.directGainExport = directGainExport;
    result.h0_modal = h0_modal;
    result.h0_display = h0_display;
    result.slope0_per_second = slope0_per_second;
    result.slope0_per_sample = slope0_per_sample;
    result.slope0_error_from_minus_one = slope0_error_from_minus_one;
    result.exportResidueNorm = exportResidueNorm;
    result.h0_export = h0_export;
    result.slope0_export = slope0_export;
    result.slope0_export_error_from_minus_one = slope0_export_error_from_minus_one;
    result.area = area;
    result.HhpInf_before = HhpInf_before;
    result.HhpInf = HhpInf;
    result.Hhp0 = Hhp0;
end

function [b_dc, a_dc] = make_dc_blocker(dc_poles_rad)
    b_dc = 1;
    a_dc = 1;
    for i = 1:numel(dc_poles_rad)
        b_dc = conv(b_dc, [1 0]);
        a_dc = conv(a_dc, [1 dc_poles_rad(i)]);
    end
end

function dc_poles_rad = build_dc_poles(dc_block_rad, dc_block_order, dc_spread_oct)
    if isscalar(dc_block_rad)
        if dc_block_order == 1
            dc_poles_rad = dc_block_rad;
            return;
        end
        if dc_spread_oct == 0
            error('dc_spread_oct must be > 0 when scalar dc_block_rad creates multiple simple poles.');
        end

        offsets = linspace(-0.5, 0.5, dc_block_order) * dc_spread_oct;
        dc_poles_rad = dc_block_rad * (2 .^ offsets);
        return;
    end

    dc_poles_rad = dc_block_rad(:).';
    if numel(dc_poles_rad) ~= dc_block_order
        error('When dc_block_rad is a vector, its length must match dc_block_order.');
    end

    sorted_poles = sort(dc_poles_rad);
    min_gap = min(abs(diff(sorted_poles)));
    if min_gap <= 1e-9 * max(1, max(abs(sorted_poles)))
        error('dc_block_rad contains repeated or near-repeated poles. Use distinct simple poles for OnePoleModal.');
    end
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

function [b, a] = cancel_origin_common_factors(b, a)
    b = b(:).';
    a = a(:).';
    scale = max([1, abs(b), abs(a)]);
    tol = 1e-12 * scale;

    while numel(b) > 1 && numel(a) > 1 && abs(b(end)) <= tol && abs(a(end)) <= tol
        b = b(1:end-1);
        a = a(1:end-1);
    end

    if abs(a(end)) <= tol
        warning('Denominator still has an uncancelled pole at the origin.');
    end
end

function v = trim_leading_zeros(v)
    v = v(:).';
    % Analog polynomials can span enormous coefficient ranges. Do not use a
    % relative-to-max tolerance here, or the valid leading coefficient 1 can
    % be mistaken for zero.
    idx = find(v ~= 0, 1, 'first');
    if isempty(idx)
        v = 0;
    else
        v = v(idx:end);
    end
end

function h0 = initial_value(b, a)
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);
    deg_b = numel(b) - 1;
    deg_a = numel(a) - 1;

    if deg_a == deg_b + 1
        h0 = b(1) / a(1);
    elseif deg_a > deg_b + 1
        h0 = 0;
    else
        h0 = Inf;
    end
end

function slope0 = initial_derivative(b, a)
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);
    deg_b = numel(b) - 1;
    deg_a = numel(a) - 1;

    if deg_a == deg_b + 2
        slope0 = b(1) / a(1);
    elseif deg_a > deg_b + 2
        slope0 = 0;
    else
        slope0 = Inf;
    end
end

function gain = high_frequency_gain(b, a)
    b = trim_leading_zeros(b);
    a = trim_leading_zeros(a);
    deg_b = numel(b) - 1;
    deg_a = numel(a) - 1;

    if deg_b == deg_a
        gain = b(1) / a(1);
    elseif deg_b < deg_a
        gain = 0;
    else
        gain = Inf;
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
