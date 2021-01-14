%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3.1

% Part a
mac_rgb_image = imread('assets/maccropped.jpg');
mac_gs_image = rgb2gray(mac_rgb_image);

imshow(mac_gs_image);
%%
% Part b

sobel_vertical = [
            -1 0 1; 
            -2 0 2; 
            -1 0 1;
          ];

sobel_horizontal = [
           -1 -2 -1; 
            0  0  0; 
            1  2  1;
          ];

result_vertical = conv2(mac_gs_image, sobel_vertical);
figure, imshow(uint8(result_vertical))

result_horizontal = conv2(mac_gs_image, sobel_horizontal); 
figure, imshow(uint8(result_horizontal)) 
%%
% Part c
canny_image = result_vertical.^2 + result_horizontal.^2;
result_square_rooted = sqrt(canny_image);
figure, imshow(uint8(canny_image));
figure, imshow(uint8(result_square_rooted));
%%
% Part d
threshold_list = [0, 50, 100, 150, 200, 250, 300, 350, 400];

for i = 1:9  
    thresholded_image = result_square_rooted>threshold_list(i);
    subplot(3,3,i);
    imshow(thresholded_image);
end
%%
% Part e (i)
tl = 0.04;
th = 0.1;
sigma = 1.0;

canny_image = edge(mac_gs_image, 'canny', [tl th], sigma(1));
figure, imshow(canny_image)

sigmaList = [2.0, 3.0, 4.0, 5.0];

figure
for i = 1:4  
    canny_image = edge(mac_gs_image, 'canny', [tl th], sigmaList(i));
    subplot(3,2,i);
    imshow(canny_image);
end
%%
clc
tl = [0.0001, 0.05, 0.075, 0.0999];
th = 0.1;
sigma = 1.0;

figure
for i = 1:4
    canny_image = edge(mac_gs_image, 'canny' ,[tl(i) th], sigma);
    subplot(2,2,i);
    imshow(canny_image);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3.2

% a
mac_rgb_image = imread('assets/maccropped.jpg');
mac_gs_image = rgb2gray(mac_rgb_image);
tl = 0.04; th = 0.1; sigmaList = 1.0;
canny_image = edge(mac_gs_image,'canny',[tl th],sigmaList);
imshow(canny_image);
%%
% Part b
[H, xp] = radon(canny_image);
imshow(uint8(H));

%%
% Part c
clc
imagesc(uint8(H));
colormap('default');
[x,y] = ginput(1)


% Part d
theta = round(x);
radius = xp(round(y));

[A, B] = pol2cart(theta*pi/180, radius);
B = -B;
C = A*(A+179) + B*(B+145);

% Part e
xl = 0;
yl = (C - A * xl) / B;

xr = 357;
yr = (C - A * xr) / B;

%f
imshow(mac_gs_image)
line([xl xr], [yl yr])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3.3
% Part b
left_corridor = imread('assets/corridorl.jpg'); 
left_corridor = rgb2gray(left_corridor);
figure, imshow(left_corridor);

right_corridor = imread('assets/corridorr.jpg');
right_corridor = rgb2gray(right_corridor);
figure, imshow(right_corridor);
%%
% Part c
D = stereoFunction(left_corridor, right_corridor);
figure, imshow(D,[-15 15]);

res = imread('assets/corridor_disp.jpg');
figure, imshow(res);
%%
% Part d
left_triclops = imread('assets/triclopsi2l.jpg'); 
left_triclops = rgb2gray(left_triclops);
figure, imshow(left_triclops);

right_triclops = imread('assets/triclopsi2r.jpg');
right_triclops = rgb2gray(right_triclops);
figure, imshow(right_triclops);

D = stereoFunction(left_triclops, right_triclops);
figure, imshow(D,[-15 15]);

res = imread('assets/triclopsid.jpg');
figure, imshow(res);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stereoImage = stereoFunction(image_l, image_r)

template_size_x = 11;
template_size_y = 11;

num_x = floor(template_size_x/2);
num_y = floor(template_size_y/2);
[x, y] = size(image_l);

stereoImage = ones(x - template_size_x + 1, y - template_size_y + 1);

for i = 1+num_x : x-num_x
    
    for j = 1+num_y : y-num_y
        
        current_image_right = image_l(i-num_x: i+num_x, j-num_y: j+num_y);
        current_image_left = rot90(current_image_right, 2);
        MC = j; 
        MD = inf;
 
        for k = max(1+num_y , j-14) : j
            
            T = image_r(i-num_x: i+num_x, k-num_y: k+num_y);
            current_image_right = rot90(T, 2);
            
            C1 = conv2(T, current_image_right);
            C2 = conv2(T, current_image_left);
            SSD = C1(template_size_x, template_size_y) - 2 * C2(template_size_x, template_size_y);
            
            if SSD < MD
                MD = SSD;
                MC = k;
            end
        end
        
        stereoImage(i - num_x, j - num_y) = j - MC;
    end
    
end
end