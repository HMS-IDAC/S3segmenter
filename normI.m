function J = normI(I)
    Irs=imresize(I,0.01,'nearest');
    p1 = prctile(Irs(:),0.1);
    J = I-p1;
%     Jrs= imresize(J,0.05,'nearest');
    p99 = prctile(Irs(:),99.9);
    J = J/(p99-p1);
end