function [all_f,all_n,all_lyso,check_lyso] = Detect_Lysosomes_Func(filename,step_xy,step_z)

tiff_info = imfinfo(filename); % return tiff structure, one element per image
tiff_stack = imread(filename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(filename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

all_lyso = [];
counter = [0];
%diff_array = containers.Map({},{},'UniformValues',false)
xo_previous = [];
yo_previous = [];
multi_lyso = {};
for im_idx = 1:size(tiff_stack,3)
    I = single(tiff_stack(:,:,im_idx));
    if mean(I(:)) < mean(tiff_stack(:))
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
        xo_previous = xo;
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
    i = i + 1
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

r=0.34;

check_lyso = all_lyso;
dmat = tril(pdist2(all_lyso,all_lyso,'euclidean'));
dmat(dmat == 0) = 10;
[idx,idy] = find(dmat <= 2*r);
check_lyso(idy,:) = [];

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

end