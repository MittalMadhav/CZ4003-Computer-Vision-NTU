%% 2.1 Contrast Stretching

% a. 
Pc = imread('assets/mrttrainbland.jpg');
whos Pc
P = rgb2gray(Pc);

% b. 
imshow(P);

% c.
a = min(P(:));
b = max(P(:));

% d.
P_Stretched = 255/191 * imsubtract(P,13);

figure, imshow(P_Stretched)
figure, imshow(P)

min(P_Stretched(:))
max(P_Stretched(:))

% e.
imshow(P_Stretched)

%% 2.2 Histogram Equalization

% a. 
imhist(P,10)
imhist(P,256)


% b. 
P3 = histeq(P,255);
figure, imhist(P3, 10)
figure, imhist(P3,256)
   
% c.
P3 = histeq(P,255);
figure, imhist(P3, 10)
figure, imhist(P3,256)


%% 2.3 Linear Spatial Filtering

% a. 
filter1 = fspecial('gaussian', 5, 1);
filter2 = fspecial('gaussian', 5, 2);

figure, mesh(filter1)
figure, mesh(filter2)

% b. 
NTU_GaussianNoise = imread('assets/ntugn.jpg');
imshow(NTU_GaussianNoise);

% c.
NTU_GaussianFilter1 = conv2(NTU_GaussianNoise,filter1, 'same');
NTU_GaussianFilter2 = conv2(NTU_GaussianNoise,filter2, 'same');

NTU_GaussianFilterFinal1 = uint8(NTU_GaussianFilter1);
NTU_GaussianFilterFinal2 = uint8(NTU_GaussianFilter2);

figure, imshow(NTU_GaussianFilterFinal1)
figure, imshow(NTU_GaussianFilterFinal2)

% d.
NTU_SpeckleNoise = imread('assets/ntusp.jpg');
imshow(NTU_SpeckleNoise);

% e.
NTU_SpeckleFilter1 = conv2(NTU_SpeckleNoise, filter1, 'same');
NTU_SpeckleFilter2 = conv2(NTU_SpeckleNoise, filter2, 'same');

NTU_SpeckleFilterFinal1 = uint8(NTU_SpeckleFilter1);
NTU_SpeckleFilterFinal2 = uint8(NTU_SpeckleFilter2);

figure, imshow(NTU_SpeckleFilterFinal1)
figure, imshow(NTU_SpeckleFilterFinal2)


%% 2.4 Median Filtering

%Gaussian Noise
filteredG3 = medfilt2(NTU_GaussianNoise, [3 3]);
filteredG5 = medfilt2(NTU_GaussianNoise, [5 5]);
figure, imshow(filteredG3)
figure, imshow(filteredG5)

%Speckle Noise
filteredS3 = medfilt2(NTU_SpeckleNoise, [3 3]);
filteredS5 = medfilt2(NTU_SpeckleNoise, [5 5]);
figure, imshow(filteredS3)
figure, imshow(filteredS5)

%% 2.5 Suppressing Noise Interference Patterns

% a. 
parallelNoiseImage = imread('assets/pckint.jpg');
imshow(parallelNoiseImage)

% b. 
F = fft2(parallelNoiseImage);
S = abs(F);
imagesc(fftshift(S.^0.1));
colormap('default')
figure, imagesc(fftshift(S.^0.1));
figure, imshow(parallelNoiseImage);
colormap(gray)

% c.
F = fft2(parallelNoiseImage);
S = abs(F);
figure, imagesc(S.^0.1);
[x_arr,y_arr] = ginput(2)

% d.
F = fft2(parallelNoiseImage);
F(14:22, 245:253) = 0;
F(237:245, 6:14) = 0;
S2 = abs(F);
imagesc(fftshift(S2.^0.1))
colormap('default')

% e.
filteredImage = uint8(ifft2(F));
imshow(filteredImage);

% f.

% Read in Image
primate = imread('assets/primatecaged.jpg');
primate = rgb2gray(primate);
imshow(primate)

% View image before setting peaks to zero
F = fft2(primate);
S = abs(F);
imagesc(S.^0.1)

% Get the 8 peak values using ginput(8) 
[y_arr,x_arr] = ginput(8)
y_arr = round(y_arr) % round Y values 
x_arr = round(x_arr) % round X values

% Removing noise using the x, y arrays
for i = 1:8
     F(x_arr(i)-2:x_arr(i)+2, y_arr(i)-2:y_arr(i)+2) = 0; 
end     

% View the Image in the Power Spectrum
S = abs(F);
imagesc(S.^0.1)   

% Invert the Fourier Transform Operation and display image
primateNew = ifft2(F);
primateNew = uint8(ifft2(F));
imshow(primateNew)


%% 2.6 Undoing Perspective Distortion of Planar Surface

% a. 
P = imread('assets/book.jpg');
imshow(P);

% b. 
[X Y] = ginput(4);
x = [0 210 0 210];
y = [0 0 295 295];

% c.
A = [
        [X(1) , Y(1) , 1 , 0    , 0    , 0, -x(1)*X(1), -x(1)*Y(1)];
        [0    , 0    , 0 , X(1) , Y(1) , 1, -y(1)*X(1), -y(1)*Y(1)];
        [X(2) , Y(2) , 1 , 0    , 0    , 0, -x(2)*X(2), -x(2)*Y(2)];
        [0    , 0    , 0 , X(2) , Y(2) , 1, -y(2)*X(2), -y(2)*Y(2)];
        [X(3) , Y(3) , 1 , 0    , 0    , 0, -x(3)*X(3), -x(3)*Y(3)];
        [0    , 0    , 0 , X(3) , Y(3) , 1, -y(3)*X(3), -y(3)*Y(3)];
        [X(4) , Y(4) , 1 , 0    , 0    , 0, -x(4)*X(4), -x(4)*Y(4)];
        [0    , 0    , 0 , X(4) , Y(4) , 1, -y(4)*X(4), -y(4)*Y(4)];
];
v = [
     x(1); 
     y(1); 
     x(2); 
     y(2); 
     x(3); 
     y(3); 
     x(4); 
     y(4)
    ];

u = A \ v;
U = reshape([u;1], 3, 3)'; 
w = U*[X'; Y'; ones(1,4)];
w = w ./ (ones(3,1) * w(3,:));

% d.
T = maketform('projective', U');
P2 = imtransform(P, T, 'XData', [0 210], 'YData', [0 297]);

% e.
imshow(P2);