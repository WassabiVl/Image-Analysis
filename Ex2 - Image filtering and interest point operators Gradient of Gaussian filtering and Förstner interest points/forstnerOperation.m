%to use the function
%   forstnerOperation(0.5, input_exercise2)
%Image Filtering and Interest Points

function matA = forstnerOperation(sigma, Image)

%Exercise 2 Interest points (Förstner Operation)
[Gimage,IX,IY,Img] = GoG(sigma, Image); % call the GoG filter function
[imageSx, imageSy] = size(Img);
Wn = ones(5,5); %weights of local Neighborhood N
WnR = 2; % radius of kernel
Ix2 = IX.^2;
Iy2=IY.*IY;
IxIy = IX.*IY;
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
%imshow(image);
hold off
matA = POI;
end



        
    
