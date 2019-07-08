function M = ridgelikl(I,nangs,scale)

% nangs = 16;
angle = 0:180/nangs:180-180/nangs;

[nr,nc] = size(I);
nconvs = length(angle);
D = zeros(nr,nc,nconvs);

count = 0;
for a = angle
    count = count+1;
    [mr,~] = smorlet(2,scale,a,1);
    C = conv2(I,mr,'same');
%     C = C.*(C > 0);
    D(:,:,count) = C;
end
M = normalize(max(D,[],3));

end