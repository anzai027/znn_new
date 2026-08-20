function [p,s] = score(pic,img)
p = psnr(pic,img);
s = ssim(pic,img);
end
