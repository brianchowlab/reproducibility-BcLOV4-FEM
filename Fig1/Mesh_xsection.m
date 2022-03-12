load('/Users/ikuz/Downloads/params_cell_2_10s_2.5dc.txt-nl-3D.mat')
q = pdem_C.Mesh;
model = createpde();
geometryFromMesh(model,q.Nodes,q.Elements);
m = generateMesh(model,'Hmax',1.25);
out = pdemesh(m,'FaceAlpha',0.1);
set(out.Parent.Children(3), 'Visible','off') 

patch([55,55,55,55],[0,60,60,0],[-7,-7,7,7],'black','FaceAlpha',0.35)

figure
planes.n = [1,0,0];
planes.r = [55,30,0];
TR = triangulation(q.Elements',q.Nodes');
el = freeBoundary(TR);
sec = mesh_xsections(q.Nodes',el,planes);
cyto = sec{1}{1}(:,2:3);
nuc = sec{1}{2}(:,2:3);
cyto = polyshape(cyto(:,1),cyto(:,2));
nuc = polyshape(nuc(:,1),nuc(:,2));
plot(cyto)
hold on
plot(nuc)
axis off
set(gca,'linew',2)