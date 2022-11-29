function [ b ] = f_publish_fig( f, ss )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Sep 7 2022

% the first thing is to set axis equal

if nargin == 1
  ss = 's'; % small size? tick marks are set to 10 km 
end
switch ss
case 's'
  f.Position = [0, 0, 1200, 700 ]
case 't'
  f.Position = [0, 0, 1200, 700 ]
case 'I'
  f.Position = [0, 0, 1200, 900 ]
case 'b'
  f.Position = [0, 0, 900, 800 ]
case 'r'
  f.Position = [0, 0, 1200, 600 ]
end

if length(f.Children) == 1 & strcmp(f.Children.Type,'tiledlayout')
  cax = f.Children.Children;
else 
  cax = get(f,'children');
end

for i = 1:length(cax)
  ca = cax(i)
  % convert all plots to km
  switch ca.Type
  case 'colorbar'
    ca.FontSize = 16;
    ca.Label.FontSize = 16;

  case 'axes'
    %ca.InnerPosition = [0.13 0.11 0.775 0.815 ];
    b = f_publish_axis(ca,ss);
  end
end
