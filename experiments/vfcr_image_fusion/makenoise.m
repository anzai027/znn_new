function [W,pics] = makenoise(img,kind,level,n,seed)
rng(seed,'twister')
W = zeros(numel(img),n);
pics = zeros([size(img),n],'single');
for i = 1:n
    if kind == "gaussian"
        pic = imnoise(img,'gaussian',0,level);
    elseif kind == "salt"
        pic = imnoise(img,'salt & pepper',level);
    elseif kind == "rician"
        dev = sqrt(level);
        a = img+dev*randn(size(img));
        b = dev*randn(size(img));
        pic = min(sqrt(a.^2+b.^2),1);
    else
        error('VFCR:Noise','Unknown noise type.')
    end
    W(:,i) = pic(:);
    pics(:,:,i) = single(pic);
end
end
