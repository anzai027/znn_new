clearvars
close all
clc

root = fileparts(mfilename('fullpath'));
addpath(root)
addpath(fileparts(root))
out = fullfile(root,'output');
if ~exist(out,'dir')
    mkdir(out)
end
log = fullfile(out,'run.log');
if exist(log,'file')
    delete(log)
end
diary(log)

cfg = config();
meta.time = datetime('now');
meta.matlab = version;
meta.system = computer;
meta.method = "VFCR-ZNN only";

fprintf('1/3 Running BSDS300 experiments.\n')
bsds = runbsds(root,cfg);
fprintf('2/3 Running 20/50 image experiment.\n')
count = runcount(root,cfg);
fprintf('3/3 Running BrainWeb experiments.\n')
mri = runmri(root,cfg);

diary off
save(fullfile(out,'all.mat'),'-v7.3')
fprintf('Finished: %s\n',out)
