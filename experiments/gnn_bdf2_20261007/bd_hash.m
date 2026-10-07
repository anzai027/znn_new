function code = bd_hash(path)
% 记录源码的校验值。
f = fopen(path,'rb'); assert(f>=0); guard = onCleanup(@() fclose(f));
b = fread(f,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(b); code = string(lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[])));
end
