function res = runcount(root,cfg)
data = fullfile(root,'data','SIPI');
out = fullfile(root,'output');
figs = fullfile(out,'fig');
imgs = fullfile(out,'image','count');
tabs = fullfile(out,'table');
if ~exist(figs,'dir'), mkdir(figs), end
if ~exist(imgs,'dir'), mkdir(imgs), end
if ~exist(tabs,'dir'), mkdir(tabs), end

path = getfile(data,'boat.512.tiff');
img = readimg(path);
[W,pics] = makenoise(img,"gaussian",0.05,50,cfg.seed+20);
nums = [20;50];
paperp = [20.99;30.03];
papers = [0.61;0.75];
ps = zeros(2,1);
ss = zeros(2,1);
runtime = zeros(2,1);
final = zeros(2,1);
sumx = zeros(2,1);
minx = zeros(2,1);
runs = repmat(struct,2,1);

for i = 1:2
    n = nums(i);
    A = W(:,1:n);
    G = A'*A/size(A,1);
    G = (G+G')/2;
    sol = solvevfcr(G,cfg,cfg.seed+200+i);
    fused = reshape(A*sol.x,size(img));
    [ps(i),ss(i)] = score(fused,img);
    runtime(i) = sol.run;
    final(i) = sol.final;
    sumx(i) = sol.sum;
    minx(i) = sol.min;
    runs(i).count = n;
    runs(i).G = G;
    runs(i).fused = single(fused);
    runs(i).sol = sol;
    imwrite(min(max(fused,0),1),fullfile(imgs,sprintf('vfcr_%d.png',n)))
end

T = table(nums,ps,ss,paperp,papers,runtime,final,sumx,minx, ...
    'VariableNames',{'Images','VFCRPSNR','VFCRSSIM','PaperPSNR', ...
    'PaperSSIM','Runtime','FinalError','WeightSum','MinWeight'});
writetable(T,fullfile(tabs,'count.csv'))

f = figure('Visible','on','Color','w','Position',[100,100,950,650]);
tile = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
nexttile(tile)
imshow(img)
title('Original')
nexttile(tile)
imshow(pics(:,:,1))
title('Gaussian noise, v=0.05')
nexttile(tile)
imshow(runs(1).fused)
title(sprintf('VFCR-ZNN, 20 images, %.2f dB',ps(1)))
nexttile(tile)
imshow(runs(2).fused)
title(sprintf('VFCR-ZNN, 50 images, %.2f dB',ps(2)))
savefig(f,fullfile(figs,'count.fig'))
exportgraphics(f,fullfile(imgs,'count.png'),'Resolution',180)
close(f)

f = figure('Visible','on','Color','w');
semilogy(runs(1).sol.t,max(runs(1).sol.e,eps),'LineWidth',1.3)
hold on
semilogy(runs(2).sol.t,max(runs(2).sol.e,eps),'LineWidth',1.3)
grid on
xlabel('t (s)')
ylabel('||xi(t)||_2')
title('VFCR-ZNN residuals for 20 and 50 images')
xlim([0,3])
legend({'20 images','50 images'},'Location','best')
savefig(f,fullfile(figs,'count_error.fig'))
exportgraphics(f,fullfile(imgs,'count_error.png'),'Resolution',180)
close(f)

res.source = path;
res.original = single(img);
res.noisy = pics;
res.table = T;
res.runs = runs;
end
