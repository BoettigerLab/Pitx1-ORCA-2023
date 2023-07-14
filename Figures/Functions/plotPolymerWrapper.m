%For paper figure
function [] = plotPolymerWrapper(polymer, varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'pel','boolean',true};
defaults(end+1,:) = {'neurog','boolean',false};
defaults(end+1,:) = {'save','string',''};
defaults(end+1,:) = {'name','string',''};
pars = ParseVariableArguments(varargin,defaults,mfilename);

pel = 14;   pitx = 22;  ra3 = 35;   ra4 = 46;   pen = 55;   neurog = 64;
favs = [pitx,ra3,ra4,pen];
cmap2 = [255 127 39;    %pitx
    255 210 20;         %ra3
    180 210 40;         %ra4
    33 187 70];         %pen
if pars.pel
    favs = [pel, favs];
    cmap2 = [120 120 120; cmap2];
end
if pars.neurog
    favs = [favs,neurog];
    cmap2 = [cmap2;120 60 120];
end
cmap2 = cmap2/255;

% tube color
cmap = zeros(75,3) + [120 120 120]/255;
cmap(pitx:ra3,:) = repmat([255 127 39]/255,ra3-pitx+1,1);
cmap(ra3+1:ra4,:) = repmat([255 230 80]/255,ra4-ra3,1);
cmap(ra4+1:58,:) = repmat([89 216 52]/255,58-ra4,1);
f1 = figure();
f1.Position = [0 0 1600 900];

PlotPolymerTube_TzuChiao(polymer,'maxJump',700,'method','spline',...
    'colormap',cmap,'tubeRadius',5,'sphereRadius',15,'center',false,...
    'interpPts',10, 'alpha', 1, 'lightOn', false);
PlotSpheres(polymer(favs,:),'color',cmap2,'r',30);
set(f1, 'Color', 'White');
colormap(cmap); c1 = colorbar;  axis equal;

if ~isempty(pars.save)
    saveFolder = SetFigureSavePath(pars.save);
    SaveFigure(f1,'name',pars.name,'formats',{'png','eps'},'overwrite',false);
end
end