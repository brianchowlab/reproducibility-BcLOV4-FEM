%%Cyto
folder = './params_cell_1_10s_2.5dc';

im_0s = imread([folder,'/Cell-1-10s-2.5dc0000.tif']);
ROI = ReadImageJROI('./ROIs/Cell_1_10s_2.5dc_FINAL_SELE_1.roi');
co1 = ROI.mnCoordinates;

crop = [130,440,200,720];
co1(:,1) = co1(:,1) - crop(3);
co1(:,2) = co1(:,2) - crop(1);

ROI = ReadImageJROI('./ROIs/Cell_1_10s_2.5dc_FINAL_SELE_2.roi');
co2 = ROI.mnCoordinates;
co2(:,1) = co2(:,1) - crop(3);
co2(:,2) = co2(:,2) - crop(1);

pos_roi1 = polyshape(co1);
pos_roi2 = polyshape(co2);

im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
low = 200;
high = 5245;
imagesc(im_0s,[low,high])
hold on 
plot(pos_roi1,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
plot(pos_roi2,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)


set(gca,'units','pixels') % set the axes units to pixels
x = get(gca,'position') % get the position of the axes
set(gcf,'units','pixels') % set the figure units to pixels
y = get(gcf,'position') % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap gray
%%
figure
folder = './params_cell_1_10s_10dc';
im_0s = imread([folder,'/Cell-1-10s-10dc0000.tif'])';

ROI = ReadImageJROI('./ROIs/Cell_1_10s_10dc_FINAL_SELE_1.roi');
co1 = fliplr(ROI.mnCoordinates);

ROI = ReadImageJROI('./ROIs/Cell_1_10s_10dc_FINAL_SELE_2.roi');
co2 = fliplr(ROI.mnCoordinates);

crop = [75,550,1,700];

co1(:,1) = co1(:,1) - crop(3);
co1(:,2) = co1(:,2) - crop(1);

co2(:,1) = co2(:,1) - crop(3);
co2(:,2) = co2(:,2) - crop(1);

pos_roi_1 = polyshape(co1);
pos_roi_2 = polyshape(co2);

im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));

low = 50;
high = 2400;
imagesc(im_0s,[low,high])
hold on 
plot(pos_roi_1,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
plot(pos_roi_2,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

set(gca,'units','pixels') % set the axes units to pixels
x = get(gca,'position') % get the position of the axes
set(gcf,'units','pixels') % set the figure units to pixels
y = get(gcf,'position') % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap gray
%% XZ Section 10s 2.5%
folder = './xz_params_cell_1_10s_2.5dc';
im_CF = flipud(double(imread([folder,'/201_CF_3.5.tif']))');
im_WF = flipud(double(imread([folder,'/201_WF_3.5.tif']))');
folder = './xz_params_cell_1_10s_2.5dc_no_PSF';
im_GT = flipud(double(imread([folder,'/201_CF_3.5.tif']))');
crop = [7,90,100,550];
im_CF = im_CF(crop(1):crop(2),crop(3):crop(4));
im_WF = im_WF(crop(1):crop(2),crop(3):crop(4));
im_GT = im_GT(crop(1):crop(2),crop(3):crop(4));
x = 100;
subplot(2,3,1)
imagesc(im_GT)
colormap gray
axis off
subplot(2,3,2)
imagesc(im_CF)
colormap gray
axis off
subplot(2,3,3)
imagesc(im_WF)
colormap gray
axis off
subplot(2,3,4)
plot(im_GT(:,100))
camroll(-90)
xlim([0,84])
xlim([0,84])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(2,3,5)
plot(im_CF(:,100))
camroll(-90)
xlim([0,84])
xlim([0,84])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(2,3,6)
plot(im_WF(:,100))
camroll(-90)
xlim([0,84])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%% XZ Section 10s 10%
folder = './xz_params_cell_1_10s_10dc';
im_CF = flipud(double(imread([folder,'/201_CF_3.5.tif']))');
im_WF = flipud(double(imread([folder,'/201_WF_3.5.tif']))');
folder = './xz_params_cell_1_10s_10dc_no_PSF';
im_GT = flipud(double(imread([folder,'/201_CF_3.5.tif']))');
crop = [1,110,160,700];
im_CF = im_CF(crop(1):crop(2),crop(3):crop(4));
im_WF = im_WF(crop(1):crop(2),crop(3):crop(4));
im_GT = im_GT(crop(1):crop(2),crop(3):crop(4));
x = 100;
subplot(2,3,1)
imagesc(im_GT)
colormap gray
axis off
subplot(2,3,2)
imagesc(im_CF)
colormap gray
axis off
subplot(2,3,3)
imagesc(im_WF)
colormap gray
axis off
subplot(2,3,4)
plot(im_GT(:,175))
camroll(-90)
xlim([14,102])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(2,3,5)
plot(im_CF(:,175))
camroll(-90)
xlim([14,102])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(2,3,6)
plot(im_WF(:,175))
camroll(-90)
xlim([14,102])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%%
folder = './cyto_only_params_cell_1_10s_2.5dc_no_PSF';
im_0s = flipud(fliplr(double(imread([folder,'/1.tif']))));
im_200s = flipud(fliplr(double(imread([folder,'/201.tif']))));
crop = [130,450,200,700];

im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_200s = im_200s(crop(1):crop(2),crop(3):crop(4));


subplot(2,2,1)
imshow(im_0s)
caxis([100,4372])
subplot(2,2,2)
imshow(im_200s)
caxis([100,4372])
subplot(2,2,3)
plot([0:0.1:50]-10,im_0s(193,:))
ylim([3500,4500])
set(gca,'tickdir','out')
set(gca,'linew',2)
box off
subplot(2,2,4)
plot([0:0.1:50]-10,im_200s(193,:))
ylim([3750,3900])
set(gca,'tickdir','out')
set(gca,'linew',2)
box off

%%
folder = './cyto_only_params_cell_1_10s_10dc_no_PSF';
im_0s = flipud(fliplr(double(imread([folder,'/1.tif']))))';
im_200s = flipud(fliplr(double(imread([folder,'/201.tif']))))';
crop = [75,500,50,610];
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_200s = im_200s(crop(1):crop(2),crop(3):crop(4));

subplot(2,2,1)
imshow(im_0s)
caxis([100,2652])
subplot(2,2,2)
imshow(im_200s)
caxis([100,2652])
subplot(2,2,3)
plot([0:0.1:56]-10,im_0s(225,:))
ylim([2000,3000])
set(gca,'tickdir','out')
set(gca,'linew',2)
box off
subplot(2,2,4)
plot([0:0.1:56]-10,im_200s(225,:))
ylim([2000,2150])
set(gca,'tickdir','out')
set(gca,'linew',2)
box off


%%
%Tiles
f = figure;
colormap gray


folder = './params_cell_1_10s_2.5dc';
gamma = 0.6;
im_0s = double(imread([folder,'/Cell-1-10s-2.5dc0000.tif'])).^gamma;
im_1s = double(imread([folder,'/Cell-1-10s-2.5dc0005.tif'])).^gamma;
im_3s = double(imread([folder,'/Cell-1-10s-2.5dc0025.tif'])).^gamma;
im_5s = double(imread([folder,'/Cell-1-10s-2.5dc0050.tif'])).^gamma;
im_25s = double(imread([folder,'/Cell-1-10s-2.5dc0100.tif'])).^gamma;
im_50s = double(imread([folder,'/Cell-1-10s-2.5dc0200.tif'])).^gamma;

model_0s = 1.5*double(flipud(fliplr(Tiff([folder,'/1.tif']).read()))).^gamma;
model_1s = 1.5*double(flipud(fliplr(Tiff([folder,'/6.tif']).read()))).^gamma;
model_3s = 1.5*double(flipud(fliplr(Tiff([folder,'/26.tif']).read()))).^gamma;
model_5s = 1.5*double(flipud(fliplr(Tiff([folder,'/51.tif']).read()))).^gamma;
model_25s = 1.5*double(flipud(fliplr(Tiff([folder,'/101.tif']).read()))).^gamma;
model_50s = 1.5*double(flipud(fliplr(Tiff([folder,'/201.tif']).read()))).^gamma;

ROI = ReadImageJROI('ExcitationROI_1_10s_2.5dc.roi');
co = ROI.mnCoordinates;

% Crop

crop = [130,540,200,750];
co(:,1) = co(:,1) - crop(3);
co(:,2) = co(:,2) - crop(1);

pos_roi = polyshape(co);
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_1s = im_1s(crop(1):crop(2),crop(3):crop(4));
im_3s = im_3s(crop(1):crop(2),crop(3):crop(4));
im_5s = im_5s(crop(1):crop(2),crop(3):crop(4));
im_25s = im_25s(crop(1):crop(2),crop(3):crop(4));
im_50s = im_50s(crop(1):crop(2),crop(3):crop(4));

% scale = mean(im_0s(:)) / mean(model_0s(:));
model_0s = model_0s(crop(1):crop(2),crop(3):crop(4));
% model_0s = uint8(scale * model_0s);
model_1s = model_1s(crop(1):crop(2),crop(3):crop(4));
% model_1s = uint8(scale * model_1s);
model_3s = model_3s(crop(1):crop(2),crop(3):crop(4));
% model_3s = uint8(scale * model_3s);
model_5s = model_5s(crop(1):crop(2),crop(3):crop(4));
% model_5s = uint8(scale * model_5s);
model_25s = model_25s(crop(1):crop(2),crop(3):crop(4));
% model_25s = uint8(scale * model_25s);
model_50s = model_50s(crop(1):crop(2),crop(3):crop(4));
%model_50s = uint8(scale * model_50s);

h = tiledlayout(6,6, 'Padding', 'none', 'TileSpacing', 'none',...
    'InnerPosition',[0,0,1,1],'OuterPosition',[0,0,1,1]); 

low = 200^gamma;
high =8000^gamma;
im_0s(end-15:end-8,end-67:end-8) = high;
im_1s(end-15:end-8,end-67:end-8) = high;
im_3s(end-15:end-8,end-67:end-8) = high;
im_5s(end-15:end-8,end-67:end-8) = high;
im_25s(end-15:end-8,end-67:end-8) = high;
im_50s(end-15:end-8,end-67:end-8) = high;
model_0s(end-15:end-8,end-67:end-8) = high;
model_1s(end-15:end-8,end-67:end-8) = high;
model_3s(end-15:end-8,end-67:end-8) = high;
model_5s(end-15:end-8,end-67:end-8) = high;
model_25s(end-15:end-8,end-67:end-8) = high;
model_50s(end-15:end-8,end-67:end-8) = high;
nexttile
imagesc(im_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_1s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_3s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_5s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_25s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

nexttile
imagesc(model_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_1s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_3s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_5s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_25s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

pos = f.Position;

x = size(im_0s,1);
y = size(im_0s,2);

h.OuterPosition = [0,0,1,1];
h.Position = [0,0,1,1];
%f.Position = [250,250,y*6,x*6-20];

folder = './params_cell_1_10s_10dc';
im_0s = double(imread([folder,'/Cell-1-10s-10dc0000.tif'])');
im_0s = im_0s .^ gamma;
im_1s = double(imread([folder,'/Cell-1-10s-10dc0005.tif'])');
im_1s = im_1s .^ gamma;
im_3s = double(imread([folder,'/Cell-1-10s-10dc0025.tif'])');
im_3s = im_3s .^ gamma;
im_5s = double(imread([folder,'/Cell-1-10s-10dc0050.tif'])');
im_5s = im_5s .^ gamma;
im_25s = double(imread([folder,'/Cell-1-10s-10dc0100.tif'])');
im_25s = im_25s .^ gamma;
im_50s = double(imread([folder,'/Cell-1-10s-10dc0200.tif'])');
im_50s = im_50s .^ gamma;


model_0s = 1*double(flipud(fliplr(Tiff([folder,'/1.tif']).read()))');
model_0s = model_0s .^ gamma;
model_1s = 1*double(flipud(fliplr(Tiff([folder,'/6.tif']).read()))');
model_1s = model_1s .^ gamma;
model_3s = 1*double(flipud(fliplr(Tiff([folder,'/26.tif']).read()))');
model_3s = model_3s .^ gamma;
model_5s = 1*double(flipud(fliplr(Tiff([folder,'/51.tif']).read()))');
model_5s = model_5s .^ gamma;
model_25s = 1*double(flipud(fliplr(Tiff([folder,'/101.tif']).read()))');
model_25s = model_25s .^ gamma;
model_50s = 1*double(flipud(fliplr(Tiff([folder,'/201.tif']).read()))');
model_50s = model_50s .^ gamma;

ROI = ReadImageJROI('ExcitationROI_1_10s_10dc.roi');
co = ROI.mnCoordinates;
co = fliplr(co);

% Crop

crop = [75,550,1,700];

co(:,1) = co(:,1) - crop(3);
co(:,2) = co(:,2) - crop(1);

pos_roi = polyshape(co);
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_1s = im_1s(crop(1):crop(2),crop(3):crop(4));
im_3s = im_3s(crop(1):crop(2),crop(3):crop(4));
im_5s = im_5s(crop(1):crop(2),crop(3):crop(4));
im_25s = im_25s(crop(1):crop(2),crop(3):crop(4));
im_50s = im_50s(crop(1):crop(2),crop(3):crop(4));

% scale = mean(im_0s(:)) / mean(model_0s(:));
model_0s = model_0s(crop(1):crop(2),crop(3):crop(4));
% model_0s = uint8(scale * model_0s);
model_1s = model_1s(crop(1):crop(2),crop(3):crop(4));
% model_1s = uint8(scale * model_1s);
model_3s = model_3s(crop(1):crop(2),crop(3):crop(4));
% model_3s = uint8(scale * model_3s);
model_5s = model_5s(crop(1):crop(2),crop(3):crop(4));
% model_5s = uint8(scale * model_5s);
model_25s = model_25s(crop(1):crop(2),crop(3):crop(4));
% model_25s = uint8(scale * model_25s);
model_50s = model_50s(crop(1):crop(2),crop(3):crop(4));


low = 100^gamma;
high = 3200^gamma;
im_0s(end-15:end-8,end-67:end-8) = high;
im_1s(end-15:end-8,end-67:end-8) = high;
im_3s(end-15:end-8,end-67:end-8) = high;
im_5s(end-15:end-8,end-67:end-8) = high;
im_25s(end-15:end-8,end-67:end-8) = high;
im_50s(end-15:end-8,end-67:end-8) = high;
model_0s(end-15:end-8,end-67:end-8) = high;
model_1s(end-15:end-8,end-67:end-8) = high;
model_3s(end-15:end-8,end-67:end-8) = high;
model_5s(end-15:end-8,end-67:end-8) = high;
model_25s(end-15:end-8,end-67:end-8) = high;
model_50s(end-15:end-8,end-67:end-8) = high;
nexttile
imagesc(im_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_1s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_3s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_5s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_25s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

nexttile
imagesc(model_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_1s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_3s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_5s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_25s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
%%
folder = './Line3';
im_0s = imread([folder,'/0s.tif']);
im_1s = imread([folder,'/1s.tif']);
im_3s = imread([folder,'/3s.tif']);
im_5s = imread([folder,'/5s.tif']);
im_25s = imread([folder,'/25s.tif']);
im_50s = imread([folder,'/50s.tif']);

model_0s = Tiff([folder,'/1.tif']).read();
model_1s = Tiff([folder,'/2.tif']).read();
model_3s = Tiff([folder,'/4.tif']).read();
model_5s = Tiff([folder,'/6.tif']).read();
model_25s = Tiff([folder,'/26.tif']).read();
model_50s = Tiff([folder,'/51.tif']).read();

% Crop

%crop = [35,190,40,150];%Cell 2 1s 10 dc
%crop = [20,280,35,200];%Cell 5 5s 2 dc
crop = [1,200,30,170];%Cell 1 10s 1dc
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_1s = im_1s(crop(1):crop(2),crop(3):crop(4));
im_3s = im_3s(crop(1):crop(2),crop(3):crop(4));
im_5s = im_5s(crop(1):crop(2),crop(3):crop(4));
im_25s = im_25s(crop(1):crop(2),crop(3):crop(4));
im_50s = im_50s(crop(1):crop(2),crop(3):crop(4));

% scale = mean(im_0s(:)) / mean(model_0s(:));
model_0s = model_0s(crop(1):crop(2),crop(3):crop(4));
% model_0s = uint8(scale * model_0s);
model_1s = model_1s(crop(1):crop(2),crop(3):crop(4));
% model_1s = uint8(scale * model_1s);
model_3s = model_3s(crop(1):crop(2),crop(3):crop(4));
% model_3s = uint8(scale * model_3s);
model_5s = model_5s(crop(1):crop(2),crop(3):crop(4));
% model_5s = uint8(scale * model_5s);
model_25s = model_25s(crop(1):crop(2),crop(3):crop(4));
% model_25s = uint8(scale * model_25s);
model_50s = model_50s(crop(1):crop(2),crop(3):crop(4));
%model_50s = uint8(scale * model_50s);
% 


low = 50;
high = 1800;
im_0s(end-11:end-8,end-30:end-8) = high;
im_1s(end-11:end-8,end-30:end-8) = high;
im_3s(end-11:end-8,end-30:end-8) = high;
im_5s(end-11:end-8,end-30:end-8) = high;
im_25s(end-11:end-8,end-30:end-8) = high;
im_50s(end-11:end-8,end-30:end-8) = high;
model_0s(end-11:end-8,end-30:end-8) = high;
model_1s(end-11:end-8,end-30:end-8) = high;
model_3s(end-11:end-8,end-30:end-8) = high;
model_5s(end-11:end-8,end-30:end-8) = high;
model_25s(end-11:end-8,end-30:end-8) = high;
model_50s(end-11:end-8,end-30:end-8) = high;
nexttile
imagesc(im_0s,[low,high])
nexttile
imagesc(im_1s,[low,high])
nexttile
imagesc(im_3s,[low,high])
nexttile
imagesc(im_5s,[low,high])
nexttile
imagesc(im_25s,[low,high])
nexttile
imagesc(im_50s,[low,high])

nexttile
imagesc(model_0s,[low,high])
nexttile
imagesc(model_1s,[low,high])
nexttile
imagesc(model_3s,[low,high])
nexttile
imagesc(model_5s,[low,high])
nexttile
imagesc(model_25s,[low,high])
nexttile
imagesc(model_50s,[low,high])

%%
model_50s = Tiff(['./51-2.tif']).read();
crop = [5,200,30,170];
model_50s = model_50s(crop(1):crop(2),crop(3):crop(4));
model_50s(end-11:end-8,end-30:end-8) = 3000;
a = double(model_50s);
imagesc(a,[0,3000])
imwrite(a/3000,'out.jpg');
%% Relative depletion

folder = './Raw';
im_0s = imread([folder,'/0s.tif']);
im_1s = imread([folder,'/1s.tif']);
im_3s = imread([folder,'/3s.tif']);
im_5s = imread([folder,'/5s.tif']);
im_25s = imread([folder,'/25s.tif']);
im_50s = imread([folder,'/50s.tif']);

folder = './PSF';
model_0s_p = Tiff([folder,'/1.tif']).read();
model_1s_p = Tiff([folder,'/2.tif']).read();
model_3s_p = Tiff([folder,'/4.tif']).read();
model_5s_p = Tiff([folder,'/6.tif']).read();
model_25s_p = Tiff([folder,'/26.tif']).read();
model_50s_p = Tiff([folder,'/51.tif']).read();

folder = './NoPSF';
model_0s_np = Tiff([folder,'/1.tif']).read();
model_1s_np = Tiff([folder,'/2.tif']).read();
model_3s_np = Tiff([folder,'/4.tif']).read();
model_5s_np = Tiff([folder,'/6.tif']).read();
model_25s_np = Tiff([folder,'/26.tif']).read();
model_50s_np = Tiff([folder,'/51.tif']).read();

crop = [30,175,30,170];%Cell 1 10s 1dc
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));
im_1s = im_1s(crop(1):crop(2),crop(3):crop(4));
im_3s = im_3s(crop(1):crop(2),crop(3):crop(4));
im_5s = im_5s(crop(1):crop(2),crop(3):crop(4));
im_25s = im_25s(crop(1):crop(2),crop(3):crop(4));
im_50s = im_50s(crop(1):crop(2),crop(3):crop(4));

model_0s_p = model_0s_p(crop(1):crop(2),crop(3):crop(4));
model_1s_p = model_1s_p(crop(1):crop(2),crop(3):crop(4));
model_3s_p = model_3s_p(crop(1):crop(2),crop(3):crop(4));
model_5s_p = model_5s_p(crop(1):crop(2),crop(3):crop(4));
model_25s_p = model_25s_p(crop(1):crop(2),crop(3):crop(4));
model_50s_p = model_50s_p(crop(1):crop(2),crop(3):crop(4));

model_0s_np = model_0s_np(crop(1):crop(2),crop(3):crop(4));
model_1s_np = model_1s_np(crop(1):crop(2),crop(3):crop(4));
model_3s_np = model_3s_np(crop(1):crop(2),crop(3):crop(4));
model_5s_np = model_5s_np(crop(1):crop(2),crop(3):crop(4));
model_25s_np = model_25s_np(crop(1):crop(2),crop(3):crop(4));
model_50s_np = model_50s_np(crop(1):crop(2),crop(3):crop(4));

low = 50;
high = 2400;

im_0s(end-11:end-8,end-30:end-8) = high;
im_1s(end-11:end-8,end-30:end-8) = high;
im_3s(end-11:end-8,end-30:end-8) = high;
im_5s(end-11:end-8,end-30:end-8) = high;
im_25s(end-11:end-8,end-30:end-8) = high;
im_50s(end-11:end-8,end-30:end-8) = high;
model_0s_p(end-11:end-8,end-30:end-8) = high*10;
model_1s_p(end-11:end-8,end-30:end-8) = high*10;
model_3s_p(end-11:end-8,end-30:end-8) = high*10;
model_5s_p(end-11:end-8,end-30:end-8) = high*10;
model_25s_p(end-11:end-8,end-30:end-8) = high*10;
model_50s_p(end-11:end-8,end-30:end-8) = high*10;
model_0s_np(end-11:end-8,end-30:end-8) = high*10;
model_1s_np(end-11:end-8,end-30:end-8) = high*10;
model_3s_np(end-11:end-8,end-30:end-8) = high*10;
model_5s_np(end-11:end-8,end-30:end-8) = high*10;
model_25s_np(end-11:end-8,end-30:end-8) = high*10;
model_50s_np(end-11:end-8,end-30:end-8) = high*10;

f=figure
h = tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'none',...
    'InnerPosition',[0,0,1,1],'OuterPosition',[0,0,1,1]); 



nexttile
imagesc(im_50s,[low,high])
nexttile
imagesc(im_50s,[750,2000])

nexttile
imagesc(model_50s_np,[low,high])

%imagesc(log10(double(model_0s_np)),[3,4])
nexttile
imagesc(model_50s_np,[3450,3700])
%imagesc(log10(double(model_50s_np)),[3,4])
nexttile
imagesc(model_50s_p,[low,high])
nexttile
imagesc(model_50s_p,[750,2000])

pos = f.Position;

x = size(im_0s,1);
y = size(im_0s,2);

h.OuterPosition = [0,0,1,1];
h.Position = [0,0,1,1];
f.Position = [250,250,y*2,x*3-20];

%% Membrane/Cyto images

%Tiles
f = figure;
colormap jet


folder = './';
gamma = 0.6;

im_0s = 1*double(flipud(fliplr(Tiff([folder,'/c_1.tif']).read()))).^gamma;
im_50s = 1*double(flipud(fliplr(Tiff([folder,'/c_201.tif']).read()))).^gamma;


model_0s = 1*double(flipud(fliplr(Tiff([folder,'/m_1.tif']).read()))).^gamma;
model_50s = 1*double(flipud(fliplr(Tiff([folder,'/m_201.tif']).read()))).^gamma;

ROI = ReadImageJROI('ExcitationROI_1_10s_2.5dc.roi');
co = ROI.mnCoordinates;

% Crop

crop = [130,540,200,750];
co(:,1) = co(:,1) - crop(3);
co(:,2) = co(:,2) - crop(1);

pos_roi = polyshape(co);
im_0s = im_0s(crop(1):crop(2),crop(3):crop(4));

im_50s = im_50s(crop(1):crop(2),crop(3):crop(4));

model_0s = model_0s(crop(1):crop(2),crop(3):crop(4));

model_50s = model_50s(crop(1):crop(2),crop(3):crop(4));


h = tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'none',...
    'InnerPosition',[0,0,1,1],'OuterPosition',[0,0,1,1]); 

low = (200^gamma)/4;
high =(6000^gamma)/4;
im_0s(end-15:end-8,end-67:end-8) = high;

im_50s(end-15:end-8,end-67:end-8) = high;

model_0s(end-15:end-8,end-67:end-8) = high;

model_50s(end-15:end-8,end-67:end-8) = high;

nexttile
imagesc(im_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(im_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
colorbar


nexttile
imagesc(model_0s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
nexttile
imagesc(model_50s,[low,high])
axis off
hold on 
plot(pos_roi,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)
colorbar
pos = f.Position;

%x = size(im_0s,1);
%y = size(im_0s,2);

%h.OuterPosition = [0,0,1,1];
%h.Position = [0,0,1,1];
