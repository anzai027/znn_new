function res = runmri(root,cfg)
data = fullfile(root,'data','BrainWeb');
out = fullfile(root,'output');
figs = fullfile(out,'fig');
imgs = fullfile(out,'image','mri');
tabs = fullfile(out,'table');
if ~exist(figs,'dir'), mkdir(figs), end
if ~exist(imgs,'dir'), mkdir(imgs), end
if ~exist(tabs,'dir'), mkdir(tabs), end

path = getfile(data,'t1_icbm_normal_1mm_pn0_rf0.rawb');
slice = [40;100];
kind = ["gaussian";"rician"];
level = [0.06;0.0385];
paperp = [23.74;16.24];
papers = [0.61;0.52];
ps = zeros(2,1);
ss = zeros(2,1);
noisep = zeros(2,1);
noises = zeros(2,1);
runtime = zeros(2,1);
final = zeros(2,1);
sumx = zeros(2,1);
minx = zeros(2,1);
cases = repmat(struct,2,1);

for i = 1:2
    img = readvol(path,slice(i));
    [W,pics] = makenoise(img,kind(i),level(i),50,cfg.seed+30+i);
    G = W'*W/size(W,1);
    G = (G+G')/2;
    sol = solvevfcr(G,cfg,cfg.seed+300+i);
    fused = reshape(W*sol.x,size(img));
    [noisep(i),noises(i)] = score(double(pics(:,:,1)),img);
    [ps(i),ss(i)] = score(fused,img);
    runtime(i) = sol.run;
    final(i) = sol.final;
    sumx(i) = sol.sum;
    minx(i) = sol.min;
    cases(i).slice = slice(i);
    cases(i).kind = kind(i);
    cases(i).level = level(i);
    cases(i).original = single(img);
    cases(i).noisy = pics;
    cases(i).G = G;
    cases(i).fused = single(fused);
    cases(i).sol = sol;
    imwrite(min(max(fused,0),1),fullfile(imgs,sprintf('vfcr_%d.png',slice(i))))
end

T = table(slice,kind,level,noisep,noises,ps,ss,paperp,papers, ...
    runtime,final,sumx,minx,'VariableNames',{'Slice','Type','Level', ...
    'NoisyPSNR','NoisySSIM','VFCRPSNR','VFCRSSIM','PaperPSNR', ...
    'PaperSSIM','Runtime','FinalError','WeightSum','MinWeight'});
writetable(T,fullfile(tabs,'mri.csv'))

f = figure('Visible','off','Color','w','Position',[100,100,1000,650]);
tile = tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
for i = 1:2
    nexttile(tile)
    imshow(cases(i).original)
    title(sprintf('Original slice %d',slice(i)))
    nexttile(tile)
    imshow(cases(i).noisy(:,:,1))
    title(sprintf('%s noise',kind(i)))
    nexttile(tile)
    imshow(cases(i).fused)
    title(sprintf('VFCR-ZNN, %.2f dB',ps(i)))
end
savefig(f,fullfile(figs,'mri.fig'))
exportgraphics(f,fullfile(imgs,'mri.png'),'Resolution',180)
close(f)

f = figure('Visible','off','Color','w');
semilogy(cases(1).sol.t,max(cases(1).sol.e,eps),'LineWidth',1.3)
hold on
semilogy(cases(2).sol.t,max(cases(2).sol.e,eps),'LineWidth',1.3)
grid on
xlabel('t (s)')
ylabel('||xi(t)||_2')
title('VFCR-ZNN residuals for BrainWeb')
xlim([0,3])
legend({'Slice 40','Slice 100'},'Location','best')
savefig(f,fullfile(figs,'mri_error.fig'))
exportgraphics(f,fullfile(imgs,'mri_error.png'),'Resolution',180)
close(f)

res.source = path;
res.note = "MRI noise levels are inferred because the paper does not report them.";
res.table = T;
res.cases = cases;
end
