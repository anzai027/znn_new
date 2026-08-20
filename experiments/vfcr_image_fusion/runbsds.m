function res = runbsds(root,cfg)
data = fullfile(root,'data','BSDS300');
out = fullfile(root,'output');
figs = fullfile(out,'fig');
imgs = fullfile(out,'image','bsds');
tabs = fullfile(out,'table');
if ~exist(figs,'dir'), mkdir(figs), end
if ~exist(imgs,'dir'), mkdir(imgs), end
if ~exist(tabs,'dir'), mkdir(tabs), end

path1 = getfile(data,'385028.jpg');
path2 = getfile(data,'113016.jpg');
house = readimg(path1);
horse = readimg(path2);
kind = ["gaussian";"gaussian";"gaussian";"salt";"salt";"salt"];
name = ["house";"house";"house";"horse";"horse";"horse"];
level = [0.05;0.08;0.10;0.05;0.10;0.20];
paperp = [26.84;23.90;23.00;30.68;27.25;23.51];
papers = [0.73;0.61;0.57;0.84;0.75;0.61];
n = 20;

noisep = zeros(6,1);
noises = zeros(6,1);
ps = zeros(6,1);
ss = zeros(6,1);
runtime = zeros(6,1);
final = zeros(6,1);
sumx = zeros(6,1);
minx = zeros(6,1);
cases = repmat(struct,6,1);

for i = 1:6
    if i <= 3
        img = house;
    else
        img = horse;
    end
    [W,pics] = makenoise(img,kind(i),level(i),n,cfg.seed+i);
    G = W'*W/size(W,1);
    G = (G+G')/2;
    sol = solvevfcr(G,cfg,cfg.seed+100+i);
    fused = reshape(W*sol.x,size(img));
    [noisep(i),noises(i)] = score(double(pics(:,:,1)),img);
    [ps(i),ss(i)] = score(fused,img);
    runtime(i) = sol.run;
    final(i) = sol.final;
    sumx(i) = sol.sum;
    minx(i) = sol.min;
    cases(i).name = name(i);
    cases(i).kind = kind(i);
    cases(i).level = level(i);
    cases(i).source = string(iif(i<=3,path1,path2));
    cases(i).original = single(img);
    cases(i).noisy = pics;
    cases(i).G = G;
    cases(i).fused = single(fused);
    cases(i).sol = sol;
    imwrite(min(max(fused,0),1),fullfile(imgs,sprintf('vfcr_%d.png',i)))
end

id = (1:6)';
T = table(id,name,kind,level,noisep,noises,ps,ss,paperp,papers, ...
    runtime,final,sumx,minx,'VariableNames',{'Noise','Image','Type', ...
    'Level','NoisyPSNR','NoisySSIM','VFCRPSNR','VFCRSSIM', ...
    'PaperPSNR','PaperSSIM','Runtime','FinalError','WeightSum','MinWeight'});
writetable(T,fullfile(tabs,'bsds.csv'))

drawset(cases,1:3,'BSDS300 Gaussian noise',fullfile(figs,'bsds_gaussian.fig'), ...
    fullfile(imgs,'bsds_gaussian.png'))
drawset(cases,4:6,'BSDS300 salt and pepper noise',fullfile(figs,'bsds_salt.fig'), ...
    fullfile(imgs,'bsds_salt.png'))

f = figure('Visible','on','Color','w','Position',[100,100,900,600]);
hold on
for i = 1:6
    semilogy(cases(i).sol.t,max(cases(i).sol.e,eps),'LineWidth',1.2)
end
set(gca,'YScale','log')
ylim([1e-12,10])
grid on
xlabel('t (s)')
ylabel('||xi(t)||_2')
title('VFCR-ZNN residuals for BSDS300')
xlim([0,0.6])
legend(compose('Noise %d',1:6),'Location','best')
savefig(f,fullfile(figs,'bsds_error.fig'))
exportgraphics(f,fullfile(imgs,'bsds_error.png'),'Resolution',180)
close(f)

res.table = T;
res.cases = cases;
end

function drawset(cases,ids,titletext,figpath,pngpath)
f = figure('Visible','on','Color','w','Position',[100,100,1050,720]);
tile = tiledlayout(f,3,3,'Padding','compact','TileSpacing','compact');
for j = 1:3
    i = ids(j);
    nexttile(tile)
    imshow(cases(i).original)
    title('Original')
    nexttile(tile)
    imshow(cases(i).noisy(:,:,1))
    title(sprintf('Noisy, level %.2f',cases(i).level))
    nexttile(tile)
    imshow(cases(i).fused)
    title('VFCR-ZNN')
end
title(tile,titletext)
savefig(f,figpath)
exportgraphics(f,pngpath,'Resolution',180)
close(f)
end

function value = iif(test,a,b)
if test
    value = a;
else
    value = b;
end
end
