function result = design_elliptic_modal_normalized(n, fc_hz, Rp, Rs, fs, ir_len)
% design_elliptic_modal_normalized
% 设计 n 阶模拟低通椭圆滤波器，并按你的 SystemModal / PoleModal 离散实现方式
% 自动归一化留数，使离散实现的峰值频响约为 0 dB。
%
% 默认参数:
%   n       = 8
%   fc_hz   = 24000
%   Rp      = 1
%   Rs      = 40
%   fs      = 48000
%   ir_len  = 16384
%
% 输出 result 结构体字段:
%   .poles_all_rad     连续域全部极点(rad/s)
%   .residues_all      连续域全部留数
%   .poles_pos_rad     上半平面极点(rad/s)
%   .residues_pos      上半平面留数(未归一化)
%   .residues_norm     上半平面留数(归一化后)
%   .freq_hz           频响横轴(Hz)
%   .H_before          归一化前频响
%   .H_after           归一化后频响
%   .scale             留数统一缩放因子
%   .cpp_params_rad    适配你当前 C++ 的参数 [pre, pim(rad/s), rre, rim]
%   .cpp_params_hz     如果你想改 C++ 用 Hz，可用 [pre, f_hz, rre, rim]

    if nargin < 1 || isempty(n),      n = 12; end
    if nargin < 2 || isempty(fc_hz),  fc_hz = 22000; end
    if nargin < 3 || isempty(Rp),     Rp = 1; end
    if nargin < 4 || isempty(Rs),     Rs = 120; end
    if nargin < 5 || isempty(fs),     fs = 48000; end
    if nargin < 6 || isempty(ir_len), ir_len = 16384; end

    validateattributes(n,      {'numeric'}, {'scalar','integer','positive'});
    validateattributes(fc_hz,  {'numeric'}, {'scalar','positive'});
    validateattributes(Rp,     {'numeric'}, {'scalar','positive'});
    validateattributes(Rs,     {'numeric'}, {'scalar','positive'});
    validateattributes(fs,     {'numeric'}, {'scalar','positive'});
    validateattributes(ir_len, {'numeric'}, {'scalar','integer','>',256});

    %------------------------------------------------------------
    % 1) 模拟椭圆低通原型 -> 缩放到 fc_hz
    %------------------------------------------------------------
    wc = 2*pi*fc_hz;  % rad/s
    [z, p0, k0] = ellipap(n, Rp, Rs);
    [b0, a0] = zp2tf(z, p0, k0);
    [b, a] = lp2lp(b0, a0, wc);

    [r_all, p_all, k_dir] = residue(b, a); %#ok<ASGLU>

    %------------------------------------------------------------
    % 2) 只保留上半平面极点
    %    你的 PoleModal 一个对象就代表一对共轭项
    %------------------------------------------------------------
    idx = imag(p_all) > 0;
    p_pos = p_all(idx);
    r_pos = r_all(idx);

    % 按频率从低到高排序，方便看
    [~, order] = sort(imag(p_pos));
    p_pos = p_pos(order);
    r_pos = r_pos(order);

    %------------------------------------------------------------
    % 3) 用和 C++ 一致的模态生成算法，计算离散冲激响应与频响
    %------------------------------------------------------------
    tau = 0.0;  % 对应你 BuildSpectrumDb 里 InjectImpulse(tau,1.0) 的一种情况
    h_before = modal_impulse_response_exact(p_pos, r_pos, fs, ir_len, tau);
    [freq_hz, H_before] = one_sided_fft(h_before, fs);

    mag_before = abs(H_before);
    peak_before = max(mag_before);
    if peak_before < 1e-20
        error('归一化失败：离散频响峰值过小，可能参数或实现有问题。');
    end

    %------------------------------------------------------------
    % 4) 自动归一化留数：统一乘缩放因子
    %    使该离散实现的峰值频响 = 1 (0 dB)
    %------------------------------------------------------------
    scale = 1 / peak_before;
    r_norm = r_pos * scale;

    h_after = modal_impulse_response_exact(p_pos, r_norm, fs, ir_len, tau);
    [~, H_after] = one_sided_fft(h_after, fs);
    mag_after = abs(H_after);
    peak_after = max(mag_after);

    %------------------------------------------------------------
    % 5) 打包结果
    %------------------------------------------------------------
    cpp_params_rad = zeros(numel(p_pos)*4, 1);
    cpp_params_hz  = zeros(numel(p_pos)*4, 1);

    for i = 1:numel(p_pos)
        sigma = real(p_pos(i));
        omega_rad = imag(p_pos(i));
        freq_i_hz = omega_rad / (2*pi);

        cpp_params_rad((i-1)*4 + 1) = sigma;
        cpp_params_rad((i-1)*4 + 2) = omega_rad;
        cpp_params_rad((i-1)*4 + 3) = real(r_norm(i));
        cpp_params_rad((i-1)*4 + 4) = imag(r_norm(i));

        cpp_params_hz((i-1)*4 + 1) = sigma;
        cpp_params_hz((i-1)*4 + 2) = freq_i_hz;
        cpp_params_hz((i-1)*4 + 3) = real(r_norm(i));
        cpp_params_hz((i-1)*4 + 4) = imag(r_norm(i));
    end

    result = struct();
    result.poles_all_rad  = p_all;
    result.residues_all   = r_all;
    result.poles_pos_rad  = p_pos;
    result.residues_pos   = r_pos;
    result.residues_norm  = r_norm;
    result.freq_hz        = freq_hz;
    result.H_before       = H_before;
    result.H_after        = H_after;
    result.scale          = scale;
    result.cpp_params_rad = cpp_params_rad;
    result.cpp_params_hz  = cpp_params_hz;

    %------------------------------------------------------------
    % 6) 打印
    %------------------------------------------------------------
    fprintf('====================================================\n');
    fprintf('Elliptic modal design (normalized for SystemModal)\n');
    fprintf('Order n              = %d\n', n);
    fprintf('Cutoff fc            = %.6f Hz\n', fc_hz);
    fprintf('Sample rate fs       = %.6f Hz\n', fs);
    fprintf('Passband ripple Rp   = %.6f dB\n', Rp);
    fprintf('Stopband atten Rs    = %.6f dB\n', Rs);
    fprintf('Impulse length       = %d\n', ir_len);
    fprintf('Normalization scale  = %.12g\n', scale);
    fprintf('Peak |H| before      = %.12g (%.6f dB)\n', peak_before, 20*log10(peak_before));
    fprintf('Peak |H| after       = %.12g (%.6f dB)\n', peak_after, 20*log10(peak_after));
    fprintf('====================================================\n\n');

    fprintf('上半平面极点与归一化后留数（频率用 Hz 显示）:\n');
    for i = 1:numel(p_pos)
        sigma = real(p_pos(i));
        f_hz = imag(p_pos(i)) / (2*pi);
        rr = real(r_norm(i));
        ri = imag(r_norm(i));
        fprintf(['mode(%d): sigma = %.12g, f = %.12g Hz, ', ...
                 'residue = %.12g %+.12gj\n'], ...
                 i, sigma, f_hz, rr, ri);
    end
    fprintf('\n');

    fprintf('直接可拷到当前 C++(pim 仍然用 rad/s) 的参数:\n');
    print_cpp_vector(cpp_params_rad, 'filterParams');

    %------------------------------------------------------------
    % 7) 画图
    %------------------------------------------------------------
    figure('Name','Modal response before/after normalization');
    plot(freq_hz, 20*log10(max(abs(H_before), 1e-12)), 'LineWidth', 1.2); hold on;
    plot(freq_hz, 20*log10(max(abs(H_after),  1e-12)), 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Exact discrete modal response (same logic as C++ SystemModal)');
    legend('Before normalization', 'After normalization', 'Location', 'best');
    xlim([0, fs/2]);
end

function h = modal_impulse_response_exact(poles, residues, fs, N, tau)
% 严格按你的 C++ PoleModal 逻辑生成冲激响应
%
% 对应：
%   CalcPole()
%   InjectImpulse(tau, v)
%   ProcessSample()

    Ts = 1 / fs;
    M = numel(poles);
    h = zeros(N, 1);

    % 每个 mode 的状态
    a1 = zeros(M,1);
    a2 = zeros(M,1);
    z1 = zeros(M,1);
    z2 = zeros(M,1);
    step1 = zeros(M,1);

    for i = 1:M
        p = poles(i);
        r = residues(i);

        sigma = real(p);
        omega = imag(p);

        R = exp(sigma * Ts);
        O = omega * Ts;

        a1(i) = -2 * R * cos(O);
        a2(i) = R * R;
        step1(i) = exp(p * Ts);

        dt1 = (1 - tau) * Ts;
        shift = exp(p * dt1);
        A1 = r * shift;
        A2 = A1 * step1(i);

        g1 = 2 * real(A1);
        g2 = 2 * real(A2);

        z1(i) = z1(i) + g1;
        z2(i) = z2(i) + g2 + a1(i) * g1;
    end

    for n = 1:N
        y = 0;
        for i = 1:M
            yi = z1(i);
            z1(i) = -a1(i) * yi + z2(i);
            z2(i) = -a2(i) * yi;
            y = y + yi;
        end
        h(n) = y;
    end
end

function [f_hz, H] = one_sided_fft(x, fs)
    N = numel(x);
    X = fft(x);
    H = X(1:floor(N/2)+1);
    f_hz = (0:floor(N/2))' * (fs / N);
end

function print_cpp_vector(v, name)
    fprintf('std::vector<float> %s =\n{\n', name);
    for i = 1:4:numel(v)
        fprintf('    %.12gf, %.12gf, %.12gf, %.12gf', ...
            v(i), v(i+1), v(i+2), v(i+3));
        if i+3 < numel(v)
            fprintf(',');
        end
        fprintf('\n');
    end
    fprintf('};\n');
end