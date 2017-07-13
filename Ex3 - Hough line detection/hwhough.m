%to use the function
%   hwhough(0.5, input_ex3)

function matA =  hwhough(sigma, Image)
image = mat2gray(double(mean(Image, 3))); %create a double type matrix from greyscale image

% Step 1 create the GoG filter
r = abs(ceil(3 * sigma));
r2 = (r * 2) + 1;
sizeA = ceil(r2);
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
Cy = Cx.';
% Create the GOG filter for Gx 
Gx = zeros(sizeA,sizeA);
for i = 1:sizeA 
    for j = 1:sizeA
        Gx(i,j) = -(Cx(i,j)*exp(-((Cx(i,j)^2 + (Cy(i,j))^2)/(2*sigma^2))))/(2*pi*sigma^4);
    end    
end

%create the filter of Gy
Gy = transpose(Gx);
% apply the filter to the image by convolution
%using the method requied by the teacher
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

G = sqrt(Ix.^2 + Iy.^2);
%create the image from the matrix
G1= mat2gray(G);
level = graythresh(G1); %value 0.2156
% BW = edge(G,'canny');
BW = imbinarize(G1,level);

%Part2
Pmax= ceil(sqrt(imageSx.^2 + imageSy.^2));
Pind = -Pmax:Pmax;
Alpha = -90:89;
H = zeros(ceil(Pmax*2)+1,180);

% extended version
for z=1:imageSx
    for z2=1:imageSy 
        if BW(z,z2)==1
                alpha= atand(Iy(z,z2)/Ix(z,z2));
                p = z*cosd(alpha) + z2*sind(alpha);% to use the sin and cosine function in degree
                p1 = round(p);
                alpha1 = round(alpha);
                pis = find(Pind == p1);
                alphai = find(Alpha == alpha1);
                H(pis,alphai) = H(pis,alphai) + 1; 
        end
   end
end

% normal version
% for z=1:imageSx
%     for z2=1:imageSy 
%         if BW(z,z2)==1
%             for alpha = -90:89
%                 p = z*cos(alpha) + z2*sin(alpha);
%                 p = round(p);
%                 pis = find(Pind == p);
%                 alphai = find(Alpha == alpha);
%                 H(pis,alphai) = H(pis,alphai) + 1;                
%             end
%         end
%    end
% end

%Plot Hough Voting Array H
colormap(autumn(5));
imshow(H,[],'XData',Alpha,'YData',Pind,...
            'InitialMagnification','fit');
        axis on, axis normal, hold on;
P  = houghpeaks(double(H),40,'threshold',ceil(0.05*max(H(:))));
kx = Alpha(P(:,2)); ky = Pind(P(:,1));
plot(kx,ky,'s','color','white');
lines = houghlines(BW,Alpha,Pind,P,'FillGap',255,'MinLength',7);
figure, imshow(BW), hold on

for kin = 1:length(lines)
   kxy = [lines(kin).point1; lines(kin).point2];
   plot(kxy(:,1),kxy(:,2),'LineWidth',2,'Color','green'); 
   plot(kxy(1,1),kxy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(kxy(2,1),kxy(2,2),'x','LineWidth',2,'Color','red');
end


%Use Houghpeaks to find maxima of H
% maxH = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% imshow(H, [],'XData', Alpha, 'YData', Pind, 'InitialMagnification', 'fit');
% theX = Alpha(maxH(:,2));
% theY = Pind(maxH(:,1));
% axis on, axis normal, hold on;
% plot(theX, theY, 's', 'color', 'white');
% 
% %Use Houghlines to derive corresponding lines
% lines = houghlines(BW,Alpha,Pind,maxH);
% 
% figure, imshow(G1), hold on
% 
% for kx = 1:size(lines)
%     kyx = [lines(kx).point1; lines(kx).point2];
%     plot(kyx(:,1), kyx(:,2), 'LineWidth', 2, 'Color', 'green');
%     %start and end of lines
%     plot(kyx(1,1), kyx(1,2), 'x', 'LineWidth', 2, 'Color', 'red');
%     plot(kyx(2,1), kyx(2,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
% end
% hold off

matA = lines;
end		





        
    
