%to use the function
%   forstnerOperation(0.5, input_exercise2)
%Image Filtering and Interest Points

function matA = forstnerOperation(sigma, Image)
% sigma is the standard deviation

image = mat2gray(double(mean(Image, 3))); %create a double type matrix from greyscale image values between 0-1

%excercies 1 Image Filerting using (GoG) filter
% Step 1 create a filter kernel to be used
r = abs(ceil(3 * sigma)); %set kernel matrix colums\row size part 1 has to be an integer to create the matrix
r2 = (r * 2) + 1;  %set kernel matrix colums\row size part 2
sizeA = ceil(r2); % round the number up to its integer value
counter = r * - 1; 
Cx = zeros(sizeA, sizeA);
%create the Cx array
for n = 1:sizeA 
    for m = 1:sizeA
        Cx(n,m) = counter;
        counter = counter + 1;
    end
    counter = r * -1;
end
Cy = Cx.'; %Cy is the opposite direction of Cx
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

%Exercise 2 Interest points (Förstner Operation)
Wn = ones(5,5); %weights of local Neighborhood N
WnR = 2; % radius of kernel
Ix2 = Ix.^2;
Iy2=Iy.*Iy;
IxIy = Ix.*Iy;
 %convolving matrixes with Wn
%  WnIx2= conv2(Wn,Ix2);
%  WnIy2= conv2(Wn,Iy2);
%  WnIxIy= conv2(Wn,IxIy);  
 WnIx2 = zeros(imageSx,imageSy);
 WnIy2 = zeros(imageSx,imageSy);
 WnIxIy = zeros(imageSx,imageSy);
 for z = 1+WnR:imageSx-WnR
    for t = 1+WnR : imageSy-WnR
	Ix2F = Ix2(z-WnR:z+WnR,t-WnR:t+WnR);
    WnIx2(z,t) = sum(sum(Ix2F.*Wn));
    Iy2F = Iy2(z-WnR:z+WnR,t-WnR:t+WnR);
    WnIy2(z,t) = sum(sum(Iy2F.*Wn));
    IxIyF = IxIy(z-WnR:z+WnR,t-WnR:t+WnR);
    WnIxIy(z,t) = sum(sum(IxIyF.*Wn));
    end
 end
 %corneness for each point
 W = zeros(imageSx,imageSy);
 %roundness
 Q = zeros(imageSx,imageSy); 
 %Binary Mask Martrx
 Mc = zeros(imageSx,imageSy);
 %finding the points of interest
for k = 1:imageSx 
    for l = 1:imageSy
        M = [WnIx2(k,l), WnIxIy(k,l); WnIxIy(k,l), WnIy2(k,l)];
        w = ((trace(M))/2) - (sqrt((trace(M)/2)^2) - det(M)); %cornerness of the image at M
        q = (4*det(M))/((trace(M)^2).'); % Roundness of image at point M 
        if w>0
            W(k,l) = w;
            else
            W(k,l) = 0;
        end
        if q>=0 && q<=1
            Q(k,l) = q;
            else
            Q(k,l) = 0;
        end
        if q> 0.5 && q<1 && w> 0.004
            Mc(k,l) = 1;
        else
            Mc(k,l) = 0;
        end       
    end    
end
%points of interest
McW = W.*Mc;
McQ = Q.*Mc;
poi = McW.*McQ;
%poi = W.*Q;
POI = imregionalmax(poi);
[PoiR,PoiC] = find(POI);
figure
hold on
plot(PoiC,PoiR,'b+');
set(gca,'YDir','reverse');
imshow(image);
hold off
matA = POI;
end



        
    
