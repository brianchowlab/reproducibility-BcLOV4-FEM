%filename = 'Z7-cy5-1.tif';
%filename = 'Z9-cy5-2-realigned.tif';
%filename = 'Z9-cy5-2-realign.tif';
%filename = 'Z9-cy5-4.tif';
filename = 'Z9-cy5-3.tif';
tiff_info = imfinfo(filename); % return tiff structure, one element per image
tiff_stack = imread(filename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    ii
    temp_tiff = imread(filename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

step_xy = .2167;
step_z = 0.15;

%Red channel for 'Z9-cy5-1-realigned.tif' offset by +34 from IR
%Red channel for 'Z9-cy5-2-realigned.tif' offset by +20 from IR
%Red channel for 'Z9-cy5-3.tif' offset by +20 from IR
%Red channel for 'Z9-cy5-4.tif' offset by +25 from IR

%% Region max
all_lyso = [];
counter = [0];
%diff_array = containers.Map({},{},'UniformValues',false)
xo_previous = [];
yo_previous = [];
multi_lyso = {};
addpath ..
for im_idx = 1:size(tiff_stack,3)
    im_idx
    I = single(tiff_stack(:,:,im_idx));
    if mean(I(:)) < 0.5*mean(tiff_stack(:))
        continue
    end
    %I(I < (1*std(I(:))+mean(I(:)))) = 0;
    % Get the linear indices, as well as the subscripts, of the pixels
    % whose gradient magnitudes are larger than the given threshold
    [grdx, grdy] = gradient(I);
    grdmag = sqrt(grdx.^2 + grdy.^2);

    prm_grdthres=0.2*max(grdmag(:));  % change grdthres ratio to value

    grdmasklin = find(grdmag > prm_grdthres);
    [grdmask_IdxI, grdmask_IdxJ] = ind2sub(size(grdmag), grdmasklin);

    %h = imagesc(I)
    %set(h, 'AlphaData', (grdmag > prm_grdthres)+0.9);
    J = I;
    mask = grdmag > prm_grdthres;
    %se = strel('disk',1);
    %new_mask = imopen(mask,se);
    se = strel('disk',3);
    new_mask = imclose(mask,se);

    new_mask = ~new_mask;

    J(new_mask) = 0;
    %imagesc(J)
    CC = bwconncomp(J);

    x = 1:size(J,1);
    y = 1:size(J,2);
    [X,Y] = meshgrid(x,y);

    xo = [];
    yo = [];
    for i = 1:size(CC.PixelIdxList,2)
        K = zeros(size(I));
        K(CC.PixelIdxList{i}) = I(CC.PixelIdxList{i});
        [IdxI, IdxJ] = ind2sub(size(K), CC.PixelIdxList{i});
        cent_x = mean(IdxI);
        cent_y = mean(IdxJ);
        fop.StartPoint = [0,max(K(:)),cent_y,cent_x,0,max(IdxJ) - min(IdxJ),max(IdxI) - min(IdxI)];

        %imagesc(K)

        try
            [fitresult, gof] = Gauss2DRotFit(K,fop);
        catch
            continue
        end
        xc = fitresult.c1;
        yc = fitresult.c2;
        
        
        xo = [xo,xc];
        yo = [yo,yc];  

        counter_previous = counter(end):(counter(end) + size(xc,1));
        counter_previous = counter_previous(2:end);
        counter = [counter(1:end-1),counter(end):(counter(end) + size(xc,1))];


        %imagesc(fitresult(X,Y))
        %GMModel = fitgmdist(X,1);
        %BW = imregionalmax(I);
    end

    if ~isempty(xo_previous)
        d = pdist2([xo',yo'],[xo_previous',yo_previous']);
        
        [xi,yi] = find(d <= 1.01*sqrt(2));
        
        keys = [xo(xi)',yo(xi)',im_idx*ones(size(xo(xi)'))];
        vals = [xo_previous(yi)',yo_previous(yi)',(im_idx-1)*ones(size(xo_previous(yi)'))];
        
        for z=1:size(keys,1)
            multi_lyso{end+1} = [vals(z,:);keys(z,:)];
        end
        xo_previous = xo;
        yo_previous = yo;
        xo(xi) = [];
        yo(xi) = [];
        size(xi)
    else
        xo_previous = xo
        yo_previous = yo;
    end
    all_lyso = vertcat(all_lyso,[xo',yo',im_idx*ones(size(xo,2),1)]);
end


i = 2;
while i <= size(multi_lyso,2)
    for j = 1:i
        if ismember(multi_lyso{j}(end,:),multi_lyso{i}(1,:),'rows')
            multi_lyso{j} = vertcat(multi_lyso{j},multi_lyso{i}(2,:));
            multi_lyso(i) = [];
            i=2;
            break
        end
    end
    i = i + 1;
end


idx = ismember(all_lyso,cell2mat(multi_lyso'),'rows');
multi_lyso_mean = cellfun(@mean,multi_lyso,'UniformOutput',false);
all_lyso(idx,:) = [];
all_lyso = vertcat(all_lyso,cell2mat(multi_lyso_mean'));
all_lyso(:,3) = all_lyso(:,3)*step_z;
all_lyso(:,1:2) = all_lyso(:,1:2)*step_xy;

all_lyso(all_lyso(:,1)>size(tiff_stack,1)*step_xy,:) = [];
all_lyso(all_lyso(:,2)>size(tiff_stack,2)*step_xy,:) = [];
all_lyso(all_lyso(:,3)>size(tiff_stack,3)*step_z,:) = [];
all_lyso(all_lyso(:,1)<0,:) = [];
all_lyso(all_lyso(:,2)<0,:) = [];
all_lyso(all_lyso(:,3)<0,:) = [];

%plot3(all_lyso(:,1),all_lyso(:,2),all_lyso(:,3),'o')
idx_s = find(abs(all_lyso(:,3) - 19.5) < 2);
imagesc(tiff_stack(:,:,130));
hold on;
plot(all_lyso(idx_s,1)/0.2167,all_lyso(idx_s,2)/0.2167,'ro');



%%

r=0.34;

check_lyso = all_lyso;
dmat = tril(pdist2(all_lyso,all_lyso,'euclidean'));
dmat(dmat == 0) = 10;
[idx,idy] = find(dmat <= 2*r);
check_lyso(idy,:) = [];

% idx_s = find(abs(check_lyso(:,3) - 19.5) < 1);
% imagesc(tiff_stack(:,:,130));
% hold on;
% plot(check_lyso(idx_s,1)/0.2167,check_lyso(idx_s,2)/0.2167,'ro');

[x1,y1,z1] = sphere(24);
x1 = x1(:)*r;
y1 = y1(:)*r;
z1 = z1(:)*r;

figure
hold on

all_x = [];
all_y = [];
all_z = [];
all_f = [];
all_n = [];
for i = 1:size(check_lyso,1)
   x2 = x1 + check_lyso(i,1);
   y2 = y1 + check_lyso(i,2);
   z2 = z1 + check_lyso(i,3);
   
   all_x = [all_x;x2];
   all_y = [all_y;y2];
   all_z = [all_z;z2];
   
   shp = alphaShape(x2,y2,z2,1);
   f = boundaryFacets(shp);
   n = shp.Points;
   all_f = vertcat(all_f,f+size(all_n,1));
   all_n = vertcat(all_n,n);
end

t = triangulation(all_f,all_n);
trimesh(t)
%shp = alphaShape(all_x,all_y,all_z,1)
%plot(shp)
%axis equal

%%

filename = 'Z9-mCherry-4-proc.tif';
tiff_info = imfinfo(filename); % return tiff structure, one element per image
tiff_stack = imread(filename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    ii
    temp_tiff = imread(filename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

step_xy = .2167;
step_z = 0.15;

%%
mask = tiff_stack>250;
SE = strel('square',5);
SE2 = strel('square',2);
%For 'Z9-mCherry-1-proc.tif'; 300 cutoff, use 5 for SE and 3 for SE2
%For 'Z9-mCherry-4-proc.tif'; 250 cutoff, use 5 for SE and 2 for SE2
%For 'Z9-mCherry-3-proc.tif'; 230 cutoff, use 5 for SE and 2 for SE2
%For 'Z9-mCherry-2-proc.tif'; 250 cutoff, use 5 for SE and 2 for SE2



I = zeros(size(mask));
J = zeros(size(mask));
for i = 1:size(mask,3)
    J(:,:,i) = imopen(mask(:,:,i),SE2);
    J(:,:,i) = imclose(J(:,:,i),SE);
    %J(:,:,i) = imopen(mask(:,:,i),SE2);
    I(:,:,i) = imfill(J(:,:,i),'holes');
end
nuc = logical(I - J);
cyto = I;
slice=100;
green = cat(3, zeros(size(cyto(:,:,slice))),ones(size(cyto(:,:,slice))), zeros(size(cyto(:,:,slice)))); 

red = cat(3, ones(size(cyto(:,:,slice))),zeros(size(cyto(:,:,slice))), zeros(size(cyto(:,:,slice)))); 



imagesc(tiff_stack(:,:,slice))
colormap gray
hold on
h = imshow(green);
g = imshow(red);
hold off
set(h,'AlphaData', cyto(:,:,slice)*0.4)
set(g,'AlphaData', nuc(:,:,slice)*0.4)

%%

stats = regionprops3(cyto,'VoxelIdxList');
lens = cell2mat(cellfun(@(c) length(c),stats{:,'VoxelIdxList'},'UniformOutput',false));
idx= lens < 5000;
idx = cell2mat(cellfun(@(c) c,stats{idx,'VoxelIdxList'},'UniformOutput',false));
cyto_t = cyto;
cyto_t(idx) = 0;
[x,y,z] = ind2sub(size(cyto),find(cyto_t(:) == 1));


shp_c = alphaShape(x,y,z,20,'HoleThreshold',15,'RegionThreshold',1);
shp_c.RegionThreshold = 1;
[elements_c,nodes_c] = boundaryFacets(shp_c);
nodes_c(:,1) = 0.2167 * nodes_c(:,1);
nodes_c(:,2) = 0.2167 * nodes_c(:,2);
nodes_c(:,3) = 0.150 * nodes_c(:,3);
plot(shp_c)
model = createpde();
gm_c = geometryFromMesh(model,nodes_c',elements_c');
%mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',1,'Hmax',5);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');

%%
stats = regionprops3(nuc,'VoxelIdxList');
lens = cell2mat(cellfun(@(c) length(c),stats{:,'VoxelIdxList'},'UniformOutput',false));
idx= lens < 5000;
idx = cell2mat(cellfun(@(c) c,stats{idx,'VoxelIdxList'},'UniformOutput',false));
nuc_t = nuc;
nuc_t(idx) = 0;
[x,y,z] = ind2sub(size(nuc),find(nuc_t(:) == 1));


shp_n = alphaShape(x,y,z,3,'HoleThreshold',15,'RegionThreshold',1);

[elements_n,nodes_n] = boundaryFacets(shp_n);
nodes_n(:,1) = 0.2167 * nodes_n(:,1);
nodes_n(:,2) = 0.2167 * nodes_n(:,2);
nodes_n(:,3) = 0.150 * nodes_n(:,3);
plot(shp_n)

%%

nodes = [nodes_c;nodes_n];
elements = [elements_c;fliplr(elements_n)+size(nodes_c,1)];
mesh_cn = triangulation(elements,nodes);
trimesh(mesh_cn,'FaceAlpha',0.2);
%%
temp = check_lyso;
temp(:,1) = check_lyso(:,2)/0.2167;
temp(:,2) = check_lyso(:,1)/0.2167;
temp(:,3) = check_lyso(:,3)/0.15;

idx_in = inShape(shp_c,temp);

temp = temp(idx_in,:);
idx_in = inShape(shp_n,temp);
temp = temp(~idx_in,:);

plot(shp_c,'FaceAlpha',0.1)
hold on
plot3(temp(:,1),temp(:,2),temp(:,3),'ro','MarkerFaceColor', 'r')
%%

[x1,y1,z1] = sphere(8);
x1 = x1(:)*r;
y1 = y1(:)*r;
z1 = z1(:)*r;

lyso = zeros(size(temp));
lyso(:,1) = temp(:,1)*0.2167;
lyso(:,2) = temp(:,2)*0.2167;
lyso(:,3) = temp(:,3)*0.15;

figure
hold on

all_x = [];
all_y = [];
all_z = [];
all_f = [];
all_n = [];
for i = 1:size(lyso,1)
   x2 = x1 + lyso(i,1);
   y2 = y1 + lyso(i,2);
   z2 = z1 + lyso(i,3);
   
   x2t = x2 / 0.2167;
   y2t = y2 / 0.2167;
   z2t =  z2 / 0.15;
   
   idx_in = inShape(shp_c, x2t,y2t,z2t);
   idx_in_nuc = inShape(shp_n, x2t,y2t,z2t);
   
   if sum(idx_in) / length(idx_in) < 1
       continue
   end
   
   if sum(idx_in_nuc) > 0
       continue
   end
   
   all_x = [all_x;x2];
   all_y = [all_y;y2];
   all_z = [all_z;z2];
   
   shp = alphaShape(x2,y2,z2,1);
   f = boundaryFacets(shp);
   n = shp.Points;
   all_f = vertcat(all_f,f+size(all_n,1));
   all_n = vertcat(all_n,n);
end


nodes_all = [nodes;all_n];
elements_all = [elements;fliplr(all_f((0*112+1):end,:))+size(nodes,1)];
%nodes_all = [nodes];
%elements_all = [elements];

t = triangulation(elements_all,nodes_all);
trimesh(t,'FaceAlpha',0.2)
%hold on
%Ttrimesh(mesh_cn,'FaceAlpha',0.2);
%%
model = createpde();
gm_c = geometryFromMesh(model,nodes',elements');
%gm_c = geometryFromMesh(model,nodes',elements');
mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',.1,'Hmax',5);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');
pdemesh(mesh_c,'FaceAlpha',0.2)
%save([filename,'.1-5-mesh.mat'],'mesh_c','shp_n','shp_c','lyso')
%save([filename,'.1-5-mesh-nolyso.mat'],'mesh_c','shp_n','shp_c')