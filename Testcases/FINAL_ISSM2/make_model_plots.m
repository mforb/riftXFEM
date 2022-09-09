%script to make basic plots for publication
global element node
load ../FINAL_ISSM/import_issm_holly1
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
figure(1)
clf();
TR = triangulation(element,node);
cpos = TR.incenter;
