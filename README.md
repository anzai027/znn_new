# znn_new

通用VFCR-ZNN算法模块。具体的TVQP问题通过 `problem` 结构体传入，算法参数单独传入。

## 函数调用方式

函数定义为：

```matlab
function [g_dot, z_dot] = my_system( ...
    t, g, z, problem, ...
    r1, r2, lambda1, lambda2, a, p, q, delta)
```

单次调用：

```matlab
[g_dot, z_dot] = my_system( ...
    t, g, z, problem, ...
    r1, r2, lambda1, lambda2, ...
    a, p, q, delta);
```

这个函数返回当前时刻的两个导数：

```matlab
g_dot
z_dot
```

它们分别表示 \(\dot g(t)\) 和积分状态的导数 \(\dot z(t)\)，并不是完整的求解结果。

## problem 结构体

`problem` 用来保存具体TVQP问题的数据：

```matlab
problem.G  = @(t) ...;
problem.dG = @(t) ...;

problem.h  = @(t) ...;
problem.dh = @(t) ...;

problem.P  = @(t) ...;
problem.dP = @(t) ...;

problem.u  = @(t) ...;
problem.du = @(t) ...;

problem.Q  = @(t) ...;
problem.dQ = @(t) ...;

problem.v  = @(t) ...;
problem.dv = @(t) ...;
```

这里的 `@(t)` 表示该数据是时间 `t` 的函数。

如果某个矩阵不随时间变化，也要写成函数形式。例如：

```matlab
Q = -eye(2);

problem.Q  = @(t) Q;
problem.dQ = @(t) zeros(size(Q));
```

## 状态向量

设：

```matlab
t0 = 0;

n = size(problem.G(t0), 1);
m = size(problem.P(t0), 1);
w = size(problem.Q(t0), 1);
```

其中：

- `n` 是未知量 `x` 的个数；
- `m` 是等式约束的个数；
- `w` 是不等式约束的个数。

状态向量 `g` 的排列顺序必须是：

```matlab
g = [x; mu1; mu2];
```

三部分的大小分别为：

```matlab
x    % n×1
mu1  % m×1
mu2  % w×1
```

因此，`g` 和 `z` 的长度都是：

```matlab
l = n + m + w;
```

初始状态可以写成：

```matlab
x0 = zeros(n,1);
m10 = zeros(m,1);
m20 = zeros(w,1);

g0 = [x0; m10; m20];
z0 = zeros(l,1);
```

## 输入参数

| 参数 | 含义 |
|---|---|
| `t` | 当前时间 |
| `g` | 状态向量 `[x; mu1; mu2]` |
| `z` | 公式（20）中的积分状态 |
| `problem` | 具体 TVQP 问题的数据 |
| `r1`, `r2` | ZNN 模型参数 |
| `lambda1`, `lambda2` | 可变增益参数 |
| `a`, `p`, `q` | 激活函数参数 |
| `delta` | PFB 函数的扰动参数 |

## 与 ODE 求解器配合

ODE 求解器只接收一个状态向量，所以需要把 `g` 和 `z` 合并：

```matlab
Y0 = [g0; z0];

fun = @(t,Y) znn_ode( ...
    t, Y, l, problem, ...
    r1, r2, lambda1, lambda2, ...
    a, p, q, delta);

[T,Y] = ode15s(fun, tspan, Y0);
```

包装函数可以写成：

```matlab
function Y_dot = znn_ode( ...
    t, Y, l, problem, ...
    r1, r2, lambda1, lambda2, ...
    a, p, q, delta)

g = Y(1:l);
z = Y(l+1:2*l);

[g_dot, z_dot] = my_system( ...
    t, g, z, problem, ...
    r1, r2, lambda1, lambda2, ...
    a, p, q, delta);

Y_dot = [g_dot; z_dot];

end
```

求解后，`x(t)` 位于结果矩阵的前 `n` 列：

```matlab
x = Y(:,1:n);
```
