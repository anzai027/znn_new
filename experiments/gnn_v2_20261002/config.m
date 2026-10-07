function c = config()
% 保存统一的实验设置。
c.t = (0:0.005:10).';
c.delta = 1e-4;
c.l = 7;
c.rtol = 1e-5;
c.atol = 1e-7;
c.step = 0.02;
c.seconds = 180;
c.calls = 500000;
c.interval = 0.01;
c.threshold = [1e-6,1e-4];
c.dwell = 0.5;
c.tail = 5;
c.seed = 1;
c.scales = [1,10,100,1000];
c.seeds = 1:5;
c.vfcr = struct('r1',1,'r2',1,'a',1,'p',0.5,'q',0.5, ...
    'lambda1',1,'lambda2',1);
c.gnn = struct('gamma',100,'p',1);
c.pairs = [1,1;10,1;100,0.5;100,1.5];
c.primary = struct('lambda',1e-6,'eps',1e-2,'gamma',200,'rates',[100,200,200]);
end
