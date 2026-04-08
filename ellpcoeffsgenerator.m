function [r, p, k, b, a] = design_elliptic_filter(n, fc, Rp, Rs)
% design_elliptic_filter 设计 n 阶模拟低通椭圆滤波器，并输出极点与留数
%
% 输入参数（都可省略，使用默认值）:
%   n   - 滤波器阶数，默认 5
%   fc  - 截止频率(Hz)，默认 24000
%   Rp  - 通带纹波(dB)，默认 1
%   Rs  - 阻带衰减(dB)，默认 40
%
% 输出参数:
%   r   - 留数（复数）
%   p   - 极点（复数）
%   k   - 直接项
%   b,a - 传递函数分子分母系数
%
% 说明:
%   本程序设计的是模拟低通椭圆滤波器 H(s)
%   截止频率从 Hz 转换为 rad/s: wc = 2*pi*fc
%
% 示例:
%   [r,p,k,b,a] = design_elliptic_filter();
%   [r,p,k,b,a] = design_elliptic_filter(6, 24000, 0.5, 60);

    %-----------------------------
    % 1) 默认参数
    %-----------------------------
    if nargin < 1 || isempty(n)
        n = 8;
    end
    if nargin < 2 || isempty(fc)
        fc = 20000;   % Hz
    end
    if nargin < 3 || isempty(Rp)
        Rp = 1;       % dB
    end
    if nargin < 4 || isempty(Rs)
        Rs = 60;      % dB
    end

    %-----------------------------
    % 2) 参数检查
    %-----------------------------
    validateattributes(n,  {'numeric'}, {'scalar','integer','positive'}, mfilename, 'n');
    validateattributes(fc, {'numeric'}, {'scalar','positive'}, mfilename, 'fc');
    validateattributes(Rp, {'numeric'}, {'scalar','positive'}, mfilename, 'Rp');
    validateattributes(Rs, {'numeric'}, {'scalar','positive'}, mfilename, 'Rs');

    %-----------------------------
    % 3) 设计归一化模拟椭圆低通原型
    %    z, p, k 对应截止角频率 1 rad/s
    %-----------------------------
    [z, p0, k0] = ellipap(n, Rp, Rs);

    %-----------------------------
    % 4) 频率缩放到目标截止频率
    %-----------------------------
    wc = 2*pi*fc;   % rad/s

    % 将零极点增益形式转为传递函数后再进行低通频率缩放
    [b0, a0] = zp2tf(z, p0, k0);
    [b, a] = lp2lp(b0, a0, wc);

    %-----------------------------
    % 5) 求留数、极点、直接项
    %    H(s) = sum(r_i / (s - p_i)) + k
    %-----------------------------
    [r, p, k] = residue(b, a);

    %-----------------------------
    % 6) 输出结果
    %-----------------------------
    fprintf('=============================================\n');
    fprintf('模拟低通椭圆滤波器设计结果\n');
    fprintf('阶数 n              = %d\n', n);
    fprintf('截止频率 fc         = %.6f Hz\n', fc);
    fprintf('截止角频率 wc       = %.6f rad/s\n', wc);
    fprintf('通带纹波 Rp         = %.6f dB\n', Rp);
    fprintf('阻带衰减 Rs         = %.6f dB\n', Rs);
    fprintf('=============================================\n\n');

    fprintf('传递函数分子系数 b:\n');
    disp(b);

    fprintf('传递函数分母系数 a:\n');
    disp(a);

    fprintf('各极点 p_i（复数）:\n');
    for i = 1:length(p)
        fprintf('p(%d) = %.12g %+.12gj\n', i, real(p(i)), imag(p(i)));
    end
    fprintf('\n');

    fprintf('各极点对应的留数 r_i（复数）:\n');
    for i = 1:length(r)
        fprintf('r(%d) = %.12g %+.12gj\n', i, real(r(i)), imag(r(i)));
    end
    fprintf('\n');

    fprintf('直接项 k:\n');
    disp(k);

    %-----------------------------
    % 7) 可选：画幅频响应
    %-----------------------------
    figure;
    freqs(b, a);
    title(sprintf('Analog Elliptic Lowpass Filter, n=%d, fc=%.0f Hz', n, fc));
    grid on;
end