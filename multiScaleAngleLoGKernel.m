function kernel = multiScaleAngleLoGKernel(image,sigmax,sigmay,iang)
%generates an anisotropic LoG filter at angle (0, 45, 90 , 135 degrees
%selected by iang). Filter is centered and ready for fourier transform
%before convolving
[ny,nx] = size(image);

[x,y] = meshgrid(-nx/2:nx/2-1, -ny/2:ny/2-1);
x = fftshift(x);
y = fftshift(y);

nang =4;
ang = pi/nang*(iang-1);

a = (cos(ang)^2)/2/sigmax^2 + (sin(ang)^2)/2/sigmay^2;
b = -(sin(2*ang))/4/sigmax^2 + (sin(2*ang))/4/sigmay^2;
c = (sin(ang)^2)/2/sigmax^2 + (cos(ang)^2)/2/sigmay^2;
G=-exp(-a*x.^2 - 2*b*(x.*y)-c*y.^2);
d2x = ((2*a*x + 2*b*y).^2 - 2*a).*G;
d2y = ((2*b*x + 2*c*y).^2 - 2*c).*G;
kernel= sigmax*sigmay*(d2x + d2y);