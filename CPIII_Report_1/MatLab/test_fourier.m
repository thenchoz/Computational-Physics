pictest = zeros(1000,1000);
for j = 1:1000
    for k = 1:1000
        if mod(j,50) > 4 || mod(k,50) > 4
            pictest(j,k) = 1;
        end
    end
end
figure
imshow(pictest);

ftest = fftshift(fft2(pictest));
m = max(max(abs(ftest)));
ftest = ftest .* 255 ./ m;
figure
imshow(abs(ftest));

%%

pictest2 = zeros(1000,1000);
for j = 1:1000
    for k = 1:1000
        if mod(j,25) > 4 || mod(k,75) > 4
            pictest2(j,k) = 1;
        end
    end
end
figure
imshow(pictest2);

ftest2 = fftshift(fft2(pictest2));
m = max(max(abs(ftest2)));
ftest2 = ftest2 .* 255 ./ m;
figure
imshow(abs(ftest2));

%%

pictest3 = zeros(1000,1000);
for j = 1:1000
    for k = 1:1000
        if mod(j,50) > 4
            pictest3(j,k) = 1;
        end
    end
end
figure
imshow(pictest3);

ftest3 = fftshift(fft2(pictest3));
m = max(max(abs(ftest3)));
ftest3 = ftest3 .* 255 ./ m;
figure
imshow(abs(ftest3));

%%

pictest4 = zeros(1000,1000);
for j = 1:1000
    for k = 1:1000
        if mod(j,50) > 4 && mod(k,50) > 4
            pictest4(j,k) = 1;
        end
    end
end
figure
imshow(pictest4);

ftest4 = fftshift(fft2(pictest4));
m = max(max(abs(ftest4)));
ftest4 = ftest4 .* 255 ./ m;
figure
imshow(abs(ftest4));

%%

pictest5 = zeros(1000,1000);
for j = 1:1000
    for k = 1:1000
        if mod(sqrt(j^2 + k^2) ,50) > 4
            pictest5(j,k) = 1;
        end
    end
end
figure
imshow(pictest5);

ftest5 = fftshift(fft2(pictest5));
m = max(max(abs(ftest5)));
ftest5 = ftest5 .* 255 ./ m;
figure
imshow(abs(ftest5));

%%

pinit1 = imread('test.jpg');
pinit2 = imread('test2.jpg');
pinit3 = imread('test3.jpg');

pinit2 = rgb2gray(pinit2);
pinit3 = rgb2gray(pinit3);

pinit2 = double(pinit2);
pinit3 = double(pinit3);

pic1 = fftshift(fft2(pinit1));
pic2 = fftshift(fft2(pinit2));
pic3 = fftshift(fft2(pinit3));

m = max(max(abs(pic1)));
pic1 = pic1 .* 255 ./ m;
figure
imshow(pinit1);
figure
imshow(abs(pic1));
m = max(max(abs(pic2)));
pic2 = pic2 .* 255 ./ m;
figure
imshow(pinit2);
figure
imshow(abs(pic2));
m = max(max(abs(pic3)));
pic3 = pic3 .* 255 ./ m;
figure
imshow(pinit3);
figure
imshow(abs(pic3));