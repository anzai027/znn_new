function dx = ftcgnn_system(t, x, l, problem, gamma, p)
% 按论文式（9）和式（18）计算状态导数。

if ~isnumeric(gamma) || ~isscalar(gamma) || ~isreal(gamma) || ...
        ~isfinite(gamma) || gamma <= 0
    error('FTCGNN:Gamma', 'gamma 必须是正数。');
end

if ~isnumeric(p) || ~isscalar(p) || ~isreal(p) || ...
        ~isfinite(p) || p <= 0 || p >= 2
    error('FTCGNN:Exponent', 'p 必须满足 0 < p < 2。');
end

% 读取线性方程的系数。
A = problem.A(t);
b = problem.b(t);
x = x(:);
b = b(:);

if ~isnumeric(l) || ~isscalar(l) || ~isreal(l) || ...
        ~isfinite(l) || l < 1 || l ~= fix(l) || ...
        numel(x) ~= l || size(A,2) ~= l || size(A,1) ~= numel(b)
    error('FTCGNN:Size', 'x 的长度必须为 l，A 和 b 的维数必须匹配。');
end

% 计算方程残差和能量梯度。
xi = A*x-b;
s = A.'*xi;
r = norm(s,2);

% 梯度恰好为零时取零速度，避免除以零。
if r == 0
    dx = zeros(l,1);
else
    dx = -gamma*s/r^p;
end

% 式（18）的噪声直接加在状态导数上。
if isfield(problem,'noise')
    noise = problem.noise(t);
    noise = noise(:);
    if isscalar(noise)
        noise = noise*ones(l,1);
    end
    if numel(noise) ~= l
        error('FTCGNN:NoiseSize', 'noise 必须是标量或长度为 l 的向量。');
    end
    dx = dx+noise;
end

end
