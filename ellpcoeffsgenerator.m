function result = design_elliptic_modal_normalized(n, fc_hz, Rp, Rs, fs, dc_block_rad, dc_block_order, dc_spread_oct)
% Analog elliptic low-pass design + modal export for the current C++ SystemModal.
%
% Export format:
%   twoPoleParams = [pre, pim(rad/s), rre, rim, ...]
%   onePoleParams = [pre, rre, ...]
%
% DC blocker model:
%   H_dc(s) = prod_k s / (s + wd_k)
%
% Notes:
%   - n controls the elliptic low-pass order.
%   - dc_block_order controls how many DC zeros / nearby real poles are added.
%   - If dc_block_rad is a scalar and dc_block_order > 1, distinct nearby poles
%     are generated automatically around dc_block_rad so the result stays
%     compatible with simple one-pole modal sections.
%   - If dc_block_rad is a vector, its length must match dc_block_order (or, if
%     dc_block_order is omitted, the vector length becomes the order).

    if nargin < 1 || isempty(n),              n = 6; end
    if nargin < 2 || isempty(fc_hz),          fc_hz = 23000; end
    if nargin < 3 || isempty(Rp),             Rp = 100; end
    if nargin < 4 || isempty(Rs),             Rs = 1400; end
    if nargin < 5 || isempty(fs),             fs = 48000; end
    if nargin < 6 || isempty(dc_block_rad),   dc_block_rad = 100.0; end
    if nargin < 7 || isempty(dc_block_order),  dc_block_order = 0; end
    if nargin < 8 || isempty(dc_spread_oct),  dc_spread_oct = 5; end

    validateattributes(n,              {'numeric'}, {'scalar','integer','positive'});
    validateattributes(fc_hz,          {'numeric'}, {'scalar','positive'});
    validateattributes(Rp,             {'numeric'}, {'scalar','positive'});
    validateattributes(Rs,             {'numeric'}, {'scalar','positive'});
    validateattributes(fs,             {'numeric'}, {'scalar','positive'});
    validateattributes(dc_block_order, {'numeric'}, {'scalar','integer','nonnegative'});
    validateattributes(dc_spread_oct,  {'numeric'}, {'scalar','nonnegative'});

    if isscalar(dc_block_rad)
        validateattributes(dc_block_rad, {'numeric'}, {'scalar','nonnegative'});
    else
        validateattributes(dc_block_rad, {'numeric'}, {'vector','nonnegative'});
    end

    wc = 2*pi*fc_hz;

    [z, p0, k0] = ellipap(n, Rp, Rs);
    [b0, a0] = zp2tf(z, p0, k0);
    [b_lp, a_lp] = lp2lp(b0, a0, wc);

    dc_poles_rad = build_dc_poles(dc_block_rad, dc_block_order, dc_spread_oct);

    b = b_lp;
    a = a_lp;
    for i = 1:numel(dc_poles_rad)
        b = conv(b, [1 0]);
        a = conv(a, [1 dc_poles_rad(i)]);
    end

    [r_all, p_all, k_dir] = residue(b, a);

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

    if any(abs(imag(r_one)) > 1e-7 * max(1, abs(r_one)))
        warning('Some residues on numerically-real poles are not purely real. Export uses their real parts.');
    end

    H0 = polyval(b, 0) / polyval(a, 0);

    f_full_hz = unique([0, logspace(-6, log10(max(fs/2, fc_hz * 2)), 8192)]);
    H_full = freqs(b, a, 2*pi*f_full_hz);

    if isempty(dc_poles_rad)
        f_zoom_max_hz = max(fc_hz * 0.02, 1e-3);
    else
        f_zoom_max_hz = max(min(dc_poles_rad) / (2*pi) * 20, 1e-3);
    end
    f_zoom_hz = linspace(0, f_zoom_max_hz, 4096);
    H_zoom = freqs(b, a, 2*pi*f_zoom_hz);

    fprintf('====================================================\n');
    fprintf('Elliptic modal export for current SystemModal\n');
    fprintf('Elliptic order n         = %d\n', n);
    fprintf('Cutoff fc                = %.12g Hz\n', fc_hz);
    fprintf('Passband ripple Rp       = %.12g dB\n', Rp);
    fprintf('Stopband attenuation Rs  = %.12g dB\n', Rs);
    fprintf('DC blocker order         = %d\n', dc_block_order);
    if isempty(dc_poles_rad)
        fprintf('DC blocker poles         = none\n');
    else
        fprintf('DC blocker poles (rad/s) = ');
        fprintf('%.12g ', dc_poles_rad);
        fprintf('\n');
    end
    fprintf('H(0)                     = %.12g %+.12gj\n', real(H0), imag(H0));
    fprintf('|H(0)|                   = %.12g\n', abs(H0));
    fprintf('Two-pole modal count     = %d\n', numel(p_two));
    fprintf('One-pole modal count     = %d\n', numel(p_one));
    fprintf('====================================================\n\n');

    fprintf('Two-pole modes:\n');
    for i = 1:numel(p_two)
        fprintf(['two(%d): sigma = %.12g, f = %.12g Hz, ', ...
                 'residue = %.12g %+.12gj\n'], ...
                 i, real(p_two(i)), imag(p_two(i)) / (2*pi), ...
                 real(r_two(i)), imag(r_two(i)));
    end
    if isempty(p_two)
        fprintf('(none)\n');
    end
    fprintf('\n');

    fprintf('One-pole modes:\n');
    for i = 1:numel(p_one)
        fprintf('one(%d): pole = %.12g rad/s, residue = %.12g\n', ...
            i, p_one(i), real(r_one(i)));
    end
    if isempty(p_one)
        fprintf('(none)\n');
    end
    fprintf('\n');

    print_cpp_vector(twoPoleParams, 'twoPoleParams', 4);
    fprintf('\n');
    print_cpp_vector(onePoleParams, 'onePoleParams', 2);

    if any(abs(k_dir) > 1e-9)
        fprintf('\n%% WARNING: residue() returned nonzero direct term k_dir:\n');
        disp(k_dir);
    end

    figure('Name', 'Analog response debug');

    subplot(2, 1, 1);
    semilogx(f_full_hz(2:end), 20*log10(max(abs(H_full(2:end)), 1e-300)), 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Analog magnitude response');
    xlim([f_full_hz(2), f_full_hz(end)]);

    subplot(2, 1, 2);
    plot(f_zoom_hz, 20*log10(max(abs(H_zoom), 1e-300)), 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('DC zoom');
    xlim([0, f_zoom_max_hz]);

    result = struct();
    result.b = b;
    result.a = a;
    result.poles_all_rad = p_all;
    result.residues_all = r_all;
    result.k_dir = k_dir;
    result.two_poles_rad = p_two;
    result.two_residues = r_two;
    result.one_poles_rad = p_one;
    result.one_residues = r_one;
    result.twoPoleParams = twoPoleParams;
    result.onePoleParams = onePoleParams;
    result.dc_block_poles_rad = dc_poles_rad;
    result.f_full_hz = f_full_hz;
    result.H_full = H_full;
    result.f_zoom_hz = f_zoom_hz;
    result.H_zoom = H_zoom;
    result.H0 = H0;
end

function dc_poles_rad = build_dc_poles(dc_block_rad, dc_block_order, dc_spread_oct)
    if dc_block_order == 0
        dc_poles_rad = [];
        return;
    end

    if isscalar(dc_block_rad)
        if dc_block_rad <= 0
            error('dc_block_rad must be positive when dc_block_order > 0.');
        end
        if dc_block_order == 1
            dc_poles_rad = dc_block_rad;
            return;
        end
        if dc_spread_oct == 0
            error('dc_spread_oct must be > 0 when dc_block_order > 1 and dc_block_rad is scalar.');
        end

        offsets = linspace(-0.5, 0.5, dc_block_order) * dc_spread_oct;
        dc_poles_rad = dc_block_rad * (2 .^ offsets);
        return;
    end

    dc_poles_rad = dc_block_rad(:).';
    if numel(dc_poles_rad) ~= dc_block_order
        error('When dc_block_rad is a vector, its length must match dc_block_order.');
    end
    if dc_block_order > 1
        sorted_poles = sort(dc_poles_rad);
        min_gap = min(abs(diff(sorted_poles)));
        if min_gap <= 1e-9 * max(1, max(abs(sorted_poles)))
            error('dc_block_rad contains repeated or near-repeated poles that are not compatible with simple one-pole modal sections.');
        end
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
