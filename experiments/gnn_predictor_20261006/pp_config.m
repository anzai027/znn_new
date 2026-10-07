function c = pp_config()
% 保存四个模型共用的设置。
c.gamma = 200;
c.eps = 1e-2;
c.delta = 1e-4;
c.lambda = 1e-6;
c.h = 1e-5;
c.floor = 1e-10;
c.kappa = 0.05;
c.scale = [1,10,10];
c.limit = [100,5000,5000];
c.start = 0;
c.finish = 10;
c.tail = 5;
c.rtol = 1e-7;
c.atol = 1e-9;
c.step = 0.02;
c.seconds = 300;
c.calls = 500000;
c.peak = [3.4611091472136,3.4611091472136+2*pi];
c.t = (0:0.0005:10).';
for t = c.peak
    c.t = [c.t;(t-0.015:1e-5:t+0.015).';(t-0.001:1e-6:t+0.001).';t];
end
c.t = unique(c.t);
stream = RandStream('mt19937ar','Seed',1);
c.g0 = rand(stream,7,1);
end
