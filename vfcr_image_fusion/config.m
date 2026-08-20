function cfg = config()
cfg.r1 = 1;
cfg.r2 = 1;
cfg.p = 0.5;
cfg.q = 0.5;
cfg.a = 1;
cfg.l1 = 1;
cfg.l2 = 1;
cfg.delta = 1e-4;
cfg.rel = 1e-5;
cfg.abs = 1e-7;
cfg.step = 0.01;
cfg.span = linspace(0,3,1501);
cfg.tol = 1e-4;
cfg.seed = 20240813;
end
