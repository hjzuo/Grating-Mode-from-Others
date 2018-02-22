% figure;
% plot(post_offset,effectiveIndex);
figure;
eps_display=eps;
s=size(eps);
for i=1:s(1)
    for j=1:s(2)
        if eps(i,j)>=3
            eps_display(i,j)=1.484^2+(1.484^2-1.444^2);
            dss=21;
        end
    end
end
mesh(xc/1000,yc/1000,eps_display');xlabel('width/um');ylabel('height/um');axis equal
az = 0;
el = 90;
view(az, el);
