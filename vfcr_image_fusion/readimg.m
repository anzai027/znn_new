function img = readimg(path)
img = imread(path);
if size(img,3) == 3
    img = rgb2gray(img);
end
img = im2double(img);
end
