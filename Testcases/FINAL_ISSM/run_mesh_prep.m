clear element node
clear TrefineRG 
global element node
load import_issm_holly1
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
TR = triangulation(element,node);
cpos = TR.incenter;

FintX = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xx);
FintY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_yy);
FintXY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xy);
FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_H');
clear global ISSM_xx ISSM_yy ISSM_xy % without the global these only clear in this workspace!!
%figure(1)
%triplot(TR);

% WARNING: This takes time
%plotMesh(node,element,elemType,'b-','yes',figure('visible','off'))
%print([results_path,'/original_mesh'],'-dpng','-r200')

if Hidden
  f = figure('visible','off')
else
  figure()
end

indx = find(cpos(:,1)>-3.18e5);
indy = find(cpos(:,2)<-1.02e6);
in = intersect(indx,indy);
element = element(in,:);
cpos    = cpos(in,:);
f = figure();
plotMesh(node,element,elemType,'b-','yes',f)
hold on
plot(node(bc_fix,1),node(bc_fix,2),'r*')
plot(node(bc_front,1),node(bc_front,2),'cs')
print([results_path,'/section_mesh'],'-dpng','-r200')
%keyboard
clf(f);

all_node_num = unique(element);
numnode = size(all_node_num,1); 
node = node(all_node_num,:);
for i=1:3*length(element)
  element(i) = find(all_node_num==element(i));
end
edges_front2 = [];
for i=1:length(edges_front)
  if find(all_node_num==element(i,1)) & find(all_node_num==element(i,2));
    edges_front2 = [edges_front2 ; find(all_node_num==edges_front(i)),find(all_node_num==edges_front(i,2))];
  end
end
edges_front = edges_front2;
bc_fix2 = [];
for i = 1:length(bc_fix)
  if find(all_node_num==bc_fix(i))
  bc_fix2 = [ bc_fix2; find(all_node_num==bc_fix(i))];
  end
end
bc_fix = bc_fix2;

TR = triangulation(element,node);
element1 = element;
node1 = node;
cpos = TR.incenter;
f= figure();
plotMesh(node,element,elemType,'b-','yes',f)

% refinement using ameshreF
in = f_find_points_xCr(cpos,xCr,80000);
%keyboard
in1 = in;
%indx = find(cpos(:,1)>-20e3 & cpos(:,1)<110e3 );
%indy = find(cpos(:,2)>-1180e3 & cpos(:,2)<-1080e3 );

[node,element] = TrefineRG(node,element,in);

if Hidden
  figure('visible','off')
else
  figure()
end
TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement1'],'-dpng','-r200')
%keyboard
cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,25000);


%indx = find(cpos(:,1)>0e3 & cpos(:,1)<90e3 );
%indy = find(cpos(:,2)>-1160e3 & cpos(:,2)<-1100e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

if Hidden
  figure('visible','off')
else
  figure()
end
TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement2'],'-dpng','-r300')

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,20000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

if Hidden
  figure('visible','off')
else
  figure()
end
TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement3'],'-dpng','-r300');

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,5000,14000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

if Hidden
  f = figure('visible','off')
else
  f = figure()
end
TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement4'],'-dpng','-r300');
clf(); close(f);

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,3000,12000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

if Hidden
  f = figure('visible','off')
else
  f = figure()
end
TR = triangulation(element,node);
triplot(TR);
xlim([min(xs)-20000,max(xs)+20000])
ylim([min(ys)-20000,max(ys)+20000])
print([results_path,'/mesh_refinement4'],'-dpng','-r300')
clf(); close(f);

cpos = TR.incenter;
figure()
patch('faces',element,'vertices',node,'facevertexcdata',FintY(cpos)); shading flat;

%in = f_find_points_xCr(cpos,xCr,4000,6000)
%%indx = find(cpos(:,1)>39e3 & cpos(:,1)<65e3 );
%%indy = find(cpos(:,2)>-1137e3 & cpos(:,2)<-1116e3 );

%%in = intersect(indx,indy);
%[node,element] = TrefineRG(node,element,in);

all_node_num = unique(element);
numnode = size(all_node_num,1); 
node = node(all_node_num,:);
for i=1:3*length(element)
  element(i) = find(all_node_num==element(i));
end


if Hidden 
  f = figure('visible','off')
else
  f = figure()
end
numelem = size(element,1);
TR = triangulation(element,node);
triplot(TR);
xlim([min(xs)-20000,max(xs)+20000])
ylim([min(ys)-20000,max(ys)+20000])
print([results_path,'/mesh_final'],'-dpng','-r300')
clf(); close(f);

cpos = TR.incenter;
%if Hidden 
  %figure('visible','off')
%else
  %figure()
%end
%patch('faces',element,'vertices',node,'facevertexcdata',FintX(cpos)); shading flat;
%plotMesh(node,element,elemType,'b-','yes',figure())
%keyboard

bcNodes=[{bc_front} {1} {1} {bc_fix}]
edgNodes=[{edges_front} {1} {1} {1}]

