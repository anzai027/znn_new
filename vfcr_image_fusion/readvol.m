function img = readvol(path,id)
file = fopen(path,'r');
if file < 0
    error('VFCR:Data','Cannot open BrainWeb data.')
end
raw = fread(file,[181*217,181],'uint8=>double');
fclose(file);
img = reshape(raw(:,id),[181,217])'/255;
end
