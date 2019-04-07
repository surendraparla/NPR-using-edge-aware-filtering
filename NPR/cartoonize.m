function cartoon = cartoonize(image) 
sigma_s = 60;
sigma_r = 0.4;
q_level=10;
sigmaE=0.5;
phiE=1.0;
T=0.99;

shape = size(image);
image_lab = rgb2lab(image);
%image_l = 0.212*double(image(:,:,1)) + 0.715*double(image(:,:,2)) + 0.072*double(image(:,:,3));
image_l = double(image_lab(:,:,1)*(100))/255;

counter = 0;

%while (counter<2)
image_l = IC(image_l, sigma_s,sigma_r);
%counter = counter+1;
%end
I1 = quant(image_l, q_level);

while (counter<2)
image_l = NC(image_l, sigma_s,sigma_r);
counter = counter+1;
end
I2 = detectedg(image_l, sigmaE, T, phiE);
%I2 = ones(size(I1));

out_final = zeros(shape);
out_final(:,:,1) = double(I1.*I2)*255/100;
out_final(:,:,2) = image_lab(:,:,2);
out_final(:,:,3) = image_lab(:,:,3);

cartoon = lab2rgb(out_final);