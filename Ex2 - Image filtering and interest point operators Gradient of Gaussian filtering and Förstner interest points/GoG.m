%to use the function
%   GoG(0.5, input_exercise2)
%Image Filtering and Interest Points

function [Gimage,IX,IY,Img] = GoG(sigma, Image)
% sigma is the standard deviation

image = mat2gray(double(mean(Image, 3))); %create a double type matrix from greyscale image values between 0-1
Img=image;
%excercies 1 Image Filerting using (GoG) filter
% Step 1 create a filter kernel to be used
r = abs(ceil(3 * sigma)); %set kernel matrix colums\row size part 1 has to be an integer to create the matrix
r2 = (r * 2) + 1;  %set kernel matrix colums\row size part 2
sizeA = ceil(r2); % round the number up to its integer value
counter = r * - 1; %use this counter to fill the data of the matrix
Cx = zeros(sizeA, sizeA);
%create the Cx array
for n = 1:sizeA 
    for m = 1:sizeA
        Cx(n,m) = counter;
        counter = counter + 1;
    end
    counter = r * -1;
end
Cy = Cx.'; %Cy is inverse transposed of Cx
% Create the GOG filter for Gx 
Gx = zeros(sizeA,sizeA);
for i = 1:sizeA 
    for j = 1:sizeA
        Gx(i,j) = -(Cx(i,j)*exp(-((Cx(i,j)^2 + (Cy(i,j))^2)/(2*sigma^2))))/(2*pi*sigma^4);
    end    
end
%create the filter of Gy
Gy = transpose(Gx);
% apply the filter to the image by convolution (Old Method)
%Ix = conv2(Gx,image);
%Iy = conv2(Gy,image);

% apply the filter to the image by convolution (required method of convolv)
[imageSx, imageSy] = size(image);
Ix=zeros(imageSx,imageSy);
Iy=zeros(imageSx,imageSy);
for x = 1+r: imageSx-r
    for y = 1+r : imageSy-r
	f=image(x-r:x+r,y-r:y+r);
	Ix(x,y) = sum(sum(f.*Gx));
	Iy(x,y) = sum(sum(f.*Gy));
    end
end

%gradient magnitude image (g)
G = sqrt(Ix.^2 + Iy.^2);
Gimage = G;
IX = Ix;
IY=Iy;
% figure 
% imshow(mat2gray(G));
% figure
% imshow(Ix);
% figure
% imshow(Iy);
%create the image from the matrix