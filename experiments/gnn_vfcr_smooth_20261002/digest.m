function hash = digest(file)
% 计算原文件的校验值。
fid = fopen(file,'rb');
assert(fid >= 0,'无法读取文件：%s',file);
guard = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
hash = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
