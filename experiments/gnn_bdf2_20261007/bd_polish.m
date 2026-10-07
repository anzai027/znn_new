function ref = bd_polish(ref,p,delta)
% 用残差复核参考根在浮点数下的误差。
ref.shift = zeros(numel(ref.t),1); ref.xshift = ref.shift; ref.polished = ref.shift;
for k = 1:numel(ref.t)
    g = ref.g(k,:).'; old = g;
    for j = 1:3
        [xi,J] = pp_parts(ref.t(k),g,p,delta);
        step = J\xi; next = g-step;
        if norm(pp_parts(ref.t(k),next,p,delta)) > norm(xi), break, end
        g = next;
    end
    ref.shift(k) = norm(g-old); ref.xshift(k) = norm(g(1:2)-old(1:2));
    ref.polished(k) = norm(pp_parts(ref.t(k),g,p,delta));
    ref.g(k,:) = g.'; ref.xd(k,:) = g(1:2).';
end
end
