function y = vfcr_phi(x, name, a, p, q, alpha, beta)
%VFCR_PHI 计算论文列出的激活函数。

switch upper(name)
    case 'BSAF'
        c = (1 + exp(-alpha)) / (1 - exp(-alpha));
        e = exp(-alpha.*x);
        y = c .* (1 - e) ./ (1 + e);

    case 'SBPAF'
        y = 0.5 .* abs(x).^beta .* sign(x) ...
          + 0.5 .* abs(x).^(1/beta) .* sign(x);

    case 'NSBPAF'
        y = 0.5 .* x ...
          + 0.5 .* abs(x).^beta .* sign(x) ...
          + 0.5 .* abs(x).^(1/beta) .* sign(x);

    case 'NFTAF'
        y = a ...
          .* exp(abs(x).^q) ...
          .* abs(x).^p ...
          .* sign(x);

    otherwise
        error('VFCR:WrongPhi', '未知激活函数：%s。', name);
end

end
