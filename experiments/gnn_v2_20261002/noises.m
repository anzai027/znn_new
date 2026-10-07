function bank = noises(kind,seed,c)
% 提前生成随机噪声，求解过程中只读取。
bank.kind = string(kind);
bank.seed = seed;
bank.t = (0:c.interval:c.t(end)).';
bank.w = zeros(numel(bank.t)-1,c.l);
stream = RandStream('mt19937ar','Seed',seed);
if bank.kind == "gamma"
    bank.w = -log(max(rand(stream,size(bank.w)),realmin));
elseif bank.kind == "gaussian"
    bank.w = randn(stream,size(bank.w));
end
bank.random = any(bank.kind == ["gamma","gaussian"]);
end
