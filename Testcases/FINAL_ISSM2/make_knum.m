t = tiledlayout(2,2,'TileSpacing','Compact');
% tile 1
nexttile
plot([1:length(Knumerical{1,1})],Knumerical{1,1})
xlabel('step')
title('SIFs end 1')
legend({'K1','K2'})

nexttile
plot([1:length(Knumerical{1,2})],Knumerical{1,2})
xlabel('step')
title('SIFs end 2')
legend({'K1','K2'})

nexttile
plot([1:length(ThetaInc{1,1})],ThetaInc{1,1})
title('propagation angle, end 1')
xlabel('step')

nexttile
plot([1:length(ThetaInc{1,2})],ThetaInc{1,2})
title('propagation angle, end 2')
xlabel('step')
%plotMesh(node+dfa*[uxAna uyAna],element,elemType,'r-',plotNode)

figure_name = ['Knum_results'];
print([results_path,'/',figure_name],'-dpng','-r300')

