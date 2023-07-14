function PlotSpheres(lefs,varargin)
% plot lef data (3d x 2) as 3D boxes. 
% lef dimension = [#of molecules, 3 (xyz), 2(head and tail)]
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'subDivisions', 'positive', 30};
defaults(end+1,:) = {'r', 'positive', 5};
defaults(end+1,:) = {'lightingOn', 'boolean', true};
defaults(end+1,:) = {'color', 'colormap', [.3 .3 .3]};
defaults(end+1,:) = {'alpha', 'nonnegative', 1};
defaults(end+1,:) = {'method', 'string', 'wheelcart'};
% For method, choose 'left', 'right', 'dumbbell', 'wheelcart' or stick

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function


r = parameters.r;
if numel(r) == 1
    r = repmat(r,size(lefs,1),1);
end

clr = parameters.color;
if size(clr,1) == 1
    clr = repmat(clr,size(lefs,1),1);
end

alphas = parameters.alpha;
if size(alphas,1) == 1
    alphas = repmat(alphas,size(lefs,1),1);
end
    
[sx,sy,sz] = sphere(parameters.subDivisions);
meanXYZ = mean(lefs, 3);
for k = 1:size(lefs,1)
    if ~isnan(lefs(k,1))
        if strcmp(parameters.method,"left") | strcmp(parameters.method,"dumbbell")
            surf(r(k)*sx+lefs(k,1,1),r(k)*sy+lefs(k,2,1),r(k)*sz+lefs(k,3,1),...
            'EdgeColor','none','FaceColor',clr(k,:),'FaceAlpha',alphas(k)); hold on;
        end
        if strcmp(parameters.method,"right") | strcmp(parameters.method,"dumbbell")
            surf(r(k)*sx+lefs(k,1,2),r(k)*sy+lefs(k,2,2),r(k)*sz+lefs(k,3,2),...
            'EdgeColor','none','FaceColor',clr(k,:),'FaceAlpha',alphas(k)); hold on;
        end
        if strcmp(parameters.method,"wheelcart") | strcmp(parameters.method,"dumbbell") | strcmp(parameters.method,"stick")
            plot3(squeeze(lefs(k,1,:)),squeeze(lefs(k,2,:)),squeeze(lefs(k,3,:)),...
            'LineWidth', 5, 'Color',clr(k,:)); hold on;
        end
        if strcmp(parameters.method,"wheelcart")
            surf(r(k)*sx+meanXYZ(k,1),r(k)*sy+meanXYZ(k,2),r(k)*sz+meanXYZ(k,3),...
            'EdgeColor','none','FaceColor',clr(k,:),'FaceAlpha',alphas(k)); hold on;
        end 
    end
end
   
if parameters.lightingOn
    material dull;
    camlight right;
    lighting gouraud;
end