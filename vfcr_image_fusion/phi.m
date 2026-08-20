function y = phi(x,cfg)
y = cfg.a .* exp(abs(x).^cfg.q) .* abs(x).^cfg.p .* sign(x);
end
