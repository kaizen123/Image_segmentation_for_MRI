%Start script-------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Studente: Daniele Castrovilli
%%Matricola: VR386161
%%Corso di laurea magistrale: Bioinformatica e biotecnologie mediche
%%Corso: Bioimmagini AIDV
%%Analysis of multiple algorithms for segmenting an MRI image(.jpg)
%%----------
%%----------
%%1st algorithm: Simple adaptation of a segmentation algorithm on my image
%%2nd algorithm: Is made by me, a segmentation algorithm based on OTSU algorithm
%%3rd algorithm: Based on active contour modelled on my case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Polishing terminal and variables
clear
clc

%%-----------------1st algorithm
%%Using localized_seg algorithm (Quite good)
Img = imread('brain1.jpg');
m = false(size(Img));
m(170:215,200:230) = true; %Starting segmentation
seg = localized_seg(Img,m,250,1,3);%Iteratively gets near the selected feature

%%------------------2nd algorithm
%%Using OTSU algorithm for segmentation_Implemented by me
Img1 = imread('brain1.jpg');
background = imopen(Img1,strel('disk',15));figure
surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
ax = gca;
ax.YDir = 'reverse';

I2 = Img1 - background;
imshow(I2);

I3 = imadjust(I2);
imshow(I3);

bw = im2bw(I3,0.5);
figure,imshow(bw);

cc = bwconncomp(bw, 4);

%Show single component, handmade by me
grain = false(size(bw));
grain(cc.PixelIdxList{61}) = true; %61th connected object
figure,imshow(grain) %Shows the interested part of the MRI image

%Visualize connected components
labeled = labelmatrix(cc);
RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
imshow(RGB_label)
graindata = regionprops(cc,'basic'); %Area of selected obj

graindata(61).Area %Area of 61th component

%%--------------------3rd algorithm
filename = ('brain1.jpg');

img = im2double(imread(filename));
N = 100; %numero di punti del contorno

imshow(img)

[cx, cy] = ginput(1); %selecting center with user help

R=50; %radius
MAXITER = 750; %% Iteration

i = 1:N; 
theta = i*2*pi/N;
contx = cx+sin(theta)*R;
conty = cy+cos(theta)*R;

contx(N+1) = contx(1);
conty(N+1) = conty(1);
contx(N+2) = contx(2);
conty(N+2) = conty(2);

figure, imshow(img,[]); 
hold on;
plot(contx(1:N),conty(1:N),'g.','LineWidth',2); 
alpha = 3;
beta = 4;
gamma = 0;
delta = 0;
h=2*pi*R/N;

fs = fspecial('gaussian',[7 7],3);
img = imfilter(gray,fs,'replicate');

fh = fspecial('sobel');
fv = fh'; 

gy = 1*imfilter(img,fh,'replicate' );
gx = 1*imfilter(img,fv,'replicate' );
gra = sqrt(gx .* gx + gy .* gy);

imshow(gra,[0 0.5]); 
imshow(gray,[]);
hold on;

for iter = 1:MAXITER %Evolving contour
    
    %    Derivate senza cicli
    xs(2:N+1) = (contx(3:N+2)-2*contx(2:N+1)+contx(1:N))/(h*h);
    ys(2:N+1) = (conty(3:N+2)-2*conty(2:N+1)+conty(1:N))/(h*h);
    xs(1) = xs(N+1);
    ys(1) = ys(N+1);
    xs(N+2) = xs(2);
    ys(N+2) = ys(2);
    
    d(2:N+1) = sqrt((contx(3:N+2)-contx(1:N)).^2+(conty(3:N+2)-conty(1:N)).^2);
    d(1)=d(N+1);
    
    nx(2:N+1) = (conty(1:N)-conty(3:N+2))./d(2:N+1);
    ny(2:N+1) = (contx(3:N+2)-contx(1:N))./d(2:N+1);
    nx(1) = nx(N+1);
    ny(1) = ny(N+1);
    
    xss(2:N+1) = (xs(3:N+2)-2*xs(2:N+1)+xs(1:N))/(h*h);
    yss(2:N+1) = (ys(3:N+2)-2*ys(2:N+1)+ys(1:N))/(h*h);
    xss(1) = xss(N+1);
    yss(1) = yss(N+1);
    
    range =  2:N+1;
    pind =  sub2ind(round(conty(range)), round(contx(range)) );
    dx(range) = alpha * xs(range) + beta * xss(range) - gamma * gx(pind) + delta*nx(range);
    dy(range) = alpha * ys(range) + beta * yss(range)  - gamma * gy(pind) + delta*ny(range);
    
    range = 2:N+1;
    contx(range) = contx(range) + dx(range);
    conty(range) = conty(range) + dy(range);
    
    contx(1) = contx(N+1);
    conty(1) = conty(N+1);
    contx(N+2) = contx(2);
    conty(N+2) = conty(2);
end

plot(contx(1:N+1),conty(1:N+1),'r','LineWidth',1);

%%End Script--------------------------------------------------------------


