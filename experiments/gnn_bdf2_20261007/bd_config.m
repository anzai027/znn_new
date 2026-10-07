function c = bd_config()
% 保存二阶差分实验的共同条件。
c = pn_config("M8A");
c.order = 2; c.amp = 0; c.kind = "zero"; c.seed = 1;
c.predict = 1e-6;
c.rtol = 1e-11; c.atol = 1e-13;
c.seconds = 300; c.calls = 500000;
c.t = (0:0.0005:10).';
for t = c.peak
    c.t = [c.t;(t-0.015:1e-5:t+0.015).'; ...
        (t-0.001:1e-6:t+0.001).';(t-0.0001:1e-7:t+0.0001).';t];
end
c.t = unique(c.t);
c.wave = struct('w',2*pi*(20:12:200).','a',zeros(7,16),'b',zeros(7,16));
end
