function y = vfcr_phi(x, name, a, p, q, alpha, beta)
%VFCR_PHI 计算论文列出的四种激活函数。

switch upper(name) %把输入的字符串转化为大写
    case 'BSAF' %(21)
        c = (1 + exp(-alpha)) / (1 - exp(-alpha));
        e = exp(-alpha.*x);
        y = c .* (1 - e) ./ (1 + e);

    case 'SBPAF' %(22)
        y = 0.5 .* abs(x).^beta .* sign(x) ...
          + 0.5 .* abs(x).^(1/beta) .* sign(x);

    case 'NSBPAF' %(23)
        y = 0.5 .* x ...
          + 0.5 .* abs(x).^beta .* sign(x) ...
          + 0.5 .* abs(x).^(1/beta) .* sign(x);

    case 'NFTAF' %(18),Example 1传入的是这个
        y = a ...
          .* exp(abs(x).^q) ...
          .* abs(x).^p ...
          .* sign(x);

    otherwise %如果名称报错，程序会主动报错
        error('VFCR:WrongPhi', '未知激活函数：%s。', name);
end

end
