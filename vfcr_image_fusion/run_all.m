clearvars
close all
clc

root = fileparts(mfilename('fullpath'));
addpath(root)
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
meta.paper = "A Variable-Gain Fixed-Time Convergent and Robust ZNN Model for Image Fusion";

sources.bsds = "https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/";
sources.sipi = "https://sipi.usc.edu/database/";
sources.brain = "https://brainweb.bic.mni.mcgill.ca/brainweb/";

fprintf('1/3 Running BSDS300 experiments.\n')
bsds = runbsds(root,cfg);
fprintf('2/3 Running 20/50 image experiment.\n')
count = runcount(root,cfg);
fprintf('3/3 Running BrainWeb experiments.\n')
mri = runmri(root,cfg);

diary off
save(fullfile(out,'all.mat'),'-v7.3')
fprintf('Finished: %s\n',out)
