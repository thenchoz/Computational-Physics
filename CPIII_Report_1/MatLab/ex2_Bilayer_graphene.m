function ex2_Bilayer_graphene()
% ex2_Bilayer_graphene : resolve/plot second ex/report 1.
%
% Arguments: none;
%
% Returns : nothing.

% reading picture
picture = imread('bilayer.jpg');

figure
imshow(picture);

% ex2_1 : plot picture in grayscale from bilayer.jpg
% standart : from 0 to 255
picture = double(rgb2gray(picture));

figure
imshow(picture, [0, 255]);


% ex2_2 : plot picture in grayscale of fourier transform from bilayer.jpg
picture = fftshift(fft2(picture));

% translation into value from 0 to 255
m = max(max(abs(picture)));
picture = picture .* 255 ./ m;

figure
imshow(abs(picture));
% Not show between 0-255 for more contrast
    
    
% ex2_4 : combine ex_2 with color filter
filter = rgb(length(picture));

% Normally 474x474, but in case there a change of initial picture
s = size(picture);

% combine filter & picture in Fourier space
for j = 1:s(1)
    for k = 1:s(2)
        filter(j, k, :) = filter(j, k, :) .* picture(j, k);
    end
end

picture = filter;

% translation into value from 0 to 255
mR = max(max(abs(picture(:,:,1))));
mG = max(max(abs(picture(:,:,2))));
mB = max(max(abs(picture(:,:,3))));
% *5 for more contrast
picture(:,:,1) = picture(:,:,1) .* 255 .* 5 ./ mR;
picture(:,:,2) = picture(:,:,2) .* 255 .* 5 ./ mG;
picture(:,:,3) = picture(:,:,3) .* 255 .* 5 ./ mB;

figure
imshow(uint8(abs(picture)));


% ex2_5 : plot picture undo the fourier transform
% splitting into 3 layer (RGB)
pictureR = double(picture(:,:,1));
pictureG = double(picture(:,:,2));
pictureB = double(picture(:,:,3));

% undo fourier transfom
pictureR = ifft2(fftshift(pictureR));
pictureG = ifft2(fftshift(pictureG));
pictureB = ifft2(fftshift(pictureB));

% translation into value 0 to 255
mR = max(max(abs(pictureR)));
mG = max(max(abs(pictureG)));
mB = max(max(abs(pictureB)));
pictureR = pictureR .* 255 ./ mR;
pictureG = pictureG .* 255 ./ mG;
pictureB = pictureB .* 255 ./ mB;

% plot the color filter
figure
imshow(uint8(abs(picture_color(pictureR,1))), [0,255]);
figure
imshow(uint8(abs(picture_color(pictureG,2))), [0,255]);
figure
imshow(uint8(abs(picture_color(pictureB,3))), [0,255]);

% recreate the initial picture
picture(:,:,1) = pictureR;
picture(:,:,2) = pictureG;
picture(:,:,3) = pictureB;

figure
imshow(uint8(abs(picture)));


% ex2_6 : answering question
% easiest peak for each to find
[py(1), px(1)] = find(filter(:,:,1) == max(max(filter(:,:,1))));
[py(2), px(2)] = find(filter(:,:,2) == max(max(filter(:,:,2))));
[py(3), px(3)] = find(filter(:,:,3) == max(max(filter(:,:,3))));

% center of the picture
% hypotesis : centre of ring -> centre of image
c = size(filter) ./ 2;

% conversion pixel to angstrom
conv = 0.3;

r = 0;
for j = 1:3
    r = r + sqrt((px(j) - c(2))^2 + (c(1) - py(j))^2)/3;
end

lattice = conv * (c(1) * 2) / r;

disp('Q. 6');
disp(lattice);
    
    
% ex2_7 : orientation angles in Fourier space => style variable filter
% calculate angle in degree
% y = 0 at the top of the image
for j = 1:3
    angle(j) = mod(atan((c(1) - py(j))/(px(j) - c(2))) * 180 / pi,60);
end

% print results
disp('Q. 7');
disp('  Red orientation :');
disp(angle(1));
disp('  Green orientation :');
disp(angle(2));
disp('  Blue orientation :');
disp(angle(3));
disp('  Left :');
disp(angle(3) - angle(2));
disp('  Right :');
disp(angle(1) - angle(2));
end

function result = rgb(N)
% Create a filter of size NxN

% Values for the ring, depend on the purpose
c1 = N/2;
c2 = N/2;
r1 = N/8;
r2 = N/6;

% Cartesian coordinates : x,y
[x,y] = meshgrid(1:N,1:N);

% Polar coordinate : phi
phi = atan2(y-0.5*N, x-0.5*N);

% Prepare colors in hue/saturation/value (HSV) model, PI/3 => hexagon
hsv = zeros(N,N,3);
hsv(: ,: ,1) = mod(3*phi/pi ,1);

% Only a ring
for j = 1:N
    for k = 1:N
        if ((c1-j)^2 + (c2-k)^2 < r1^2) || ((c1-j)^2 + (c2-k)^2 > r2^2)
            hsv(j ,k ,2) = 0;
            hsv(j ,k ,3) = 0;
        else
            hsv(j ,k ,2) = 1;
            hsv(j ,k ,3) = 1;
        end
    end
end

% Convert to red/green/blue (RGB) model
result = hsv2rgb(hsv);
result = result * 255;

figure
imshow(uint8(result));
end

function result = picture_color(picture, RGB)
% Turn a gray scale picture into a one color (second argument) picture

result(:, :, RGB) = picture;
result(:, :, mod(RGB, 3)+ 1) = zeros(size(picture));
result(:, :, mod(RGB+1, 3)+ 1) = zeros(size(picture));
end