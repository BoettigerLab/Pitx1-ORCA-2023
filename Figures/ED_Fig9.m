% for reviewer 3 


%% Using phase-separation TADs 

%% Effect of TAD separation strengthening on E-P interaction 
folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v2_HL\'; % 
mapStep = 10;
[hl_sticky_pols,sim_maps_HL] = LoadPolySimDynamics(folder,...
           'polyStep',1,'mapStep',mapStep,'timeStep',2,'offset',0);

simMaps = cat(3,sim_maps_HL{:});
cfMap_hL = ContactFrac(simMaps,'threshold',5); 
figure(1); clf;
subplot(1,2,2); imagesc(log2(cfMap_hL)); 
axis image; colorbar;   caxis([-12,0]);
title('strong TADs')

folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v2_FL2\'; % 

[fl_sticky_pols,sim_maps_FL] = LoadPolySimDynamics(folder,...
           'polyStep',1,'mapStep',mapStep,'timeStep',2,'offset',0);

simMaps = cat(3,sim_maps_FL{:});
cfMap_fL = ContactFrac(simMaps,'threshold',5);
figure(1); 
subplot(1,2,1); imagesc(log2(cfMap_fL));
axis image; colorbar;  caxis([-12,0])
title('weak TADs')

figure(2); clf;
imagesc(cfMap_fL - cfMap_hL); axis image;
caxis([-.05 .05]);
GetColorMap('RedWhiteBlueK'); colorbar;
title('weak TADs minus strong TADs');


figure(3); clf;
imagesc(log2(cfMap_fL./cfMap_hL)); axis image;
caxis([-3 3]);
GetColorMap('RedWhiteBlueK'); colorbar;
title('log2(weak/strong)');

%%
cmap1 = [140,110,110;
        241,98,33;
        254,222,89;
        95,201,225;
        110,110,140]/255;

tads = [1,40,80,120,160,200];
nB = size(hl_sticky_pols,1);
cmap2 = zeros(nB,3);
for t=1:5
cmap2(tads(t):tads(t+1),:) = repmat(cmap1(t,:),(tads(t+1)-tads(t)+1),1);
end

figure(3); clf; 
c=1; subplot(1,2,1);
PlotPolymerTube(fl_sticky_pols{4}(1:5:end,:,1),'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2);
set(gca,'color','w'); axis off;
c=1; subplot(1,2,2);
PlotPolymerTube(hl_sticky_pols{4}(1:5:end,:,c),'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2);
set(gca,'color','w'); axis off;

figure(3); clf; 
c=3; 
subplot(1,2,1);
pol = fl_sticky_pols{4}(1:5:end,:,c);
PlotPolymerTube(pol,'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2,'center',false,'lightOn',false); hold on;
PlotSpheres(pol(tads,:),'r',1,'color','r','lightingOn',true); hold on;
set(gca,'color','w'); axis off;
 subplot(1,2,2); cla;
 c=5;
pol = hl_sticky_pols{4}(1:5:end,:,c);
PlotPolymerTube(pol,'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2,'center',false,'lightOn',false); hold on;
PlotSpheres(pol(tads,:),'r',1,'color','r','lightingOn',true); hold on;
set(gca,'color','w'); axis off;


%%
if 0  % prevent accidential overwrite execution
    f = gcf;
    name = ['U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Images\ExampleNonLE_TADs_v3.pdf'];
    exportgraphics(f,name,'ContentType','vector');
end


%% E-P bond strengthening effect on TADs
folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v2_HL\'; % 
%  folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v5_Eclustering_weak\'; 
mapStep = 10;
[hl_sticky_pols,sim_maps_HL] = LoadPolySimDynamics(folder,...
           'polyStep',1,'mapStep',mapStep,'timeStep',2,'offset',0);

simMaps = cat(3,sim_maps_HL{:});
cfMap_hL = ContactFrac(simMaps,'threshold',5); 
figure(1); clf;
subplot(1,2,1); imagesc(log2(cfMap_hL)); 
axis image; colorbar;   caxis([-12,0]);
title('weak E-P affinity');

%  folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v5_Eclustering\'; 
folder = 'T:\2023-05-05_Pitx1Sims\Model_noLE_stickyTAD_v4_Eclustering\'; %   v4 -- a bit to weak 

[fl_sticky_pols,sim_maps_FL] = LoadPolySimDynamics(folder,...
           'polyStep',1,'mapStep',mapStep,'timeStep',2,'offset',0);

simMaps = cat(3,sim_maps_FL{:});
cfMap_fL = ContactFrac(simMaps,'threshold',5);
figure(1); 
subplot(1,2,2); imagesc(log2(cfMap_fL));
title('strong E-P affinity');
axis image; colorbar;  caxis([-12,0])

figure(2); clf;
imagesc(cfMap_fL - cfMap_hL); axis image;
caxis([-.03 .03]);
colormap(flipud(GetColorMap('RedWhiteBlueK'))); colorbar;
title('weak E-P affinity minus strong E-P affinity');


figure(3); clf;
imagesc(log2(cfMap_fL./cfMap_hL)); axis image;
caxis([-3 3]);
colormap(flipud(GetColorMap('RedWhiteBlueK'))); colorbar;
title('log2(fL/hL)');

%%
cmap1 = [140,110,110;
        241,98,33;
        254,222,89;
        95,201,225;
        110,110,140]/255;

tads = [1,40,80,120,160,200];
nB = size(hl_sticky_pols,1);
cmap2 = zeros(nB,3);
for t=1:5
cmap2(tads(t):tads(t+1),:) = repmat(cmap1(t,:),(tads(t+1)-tads(t)+1),1);
end

figure(3); clf; 
c=1; subplot(1,2,1);
PlotPolymerTube(fl_sticky_pols{4}(1:5:end,:,1),'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2);
set(gca,'color','w'); axis off;
c=1; subplot(1,2,2);
PlotPolymerTube(hl_sticky_pols{4}(1:5:end,:,c),'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2);
set(gca,'color','w'); axis off;

figure(3); clf; 
c=1; 
subplot(1,2,1);
pol = fl_sticky_pols{4}(1:5:end,:,c);
PlotPolymerTube(pol,'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2,'center',false,'lightOn',false); hold on;
PlotSpheres(pol(tads,:),'r',1,'color','r','lightingOn',true); hold on;
set(gca,'color','w'); axis off;
 subplot(1,2,2); cla;
 c=6;
pol = hl_sticky_pols{4}(1:5:end,:,c);
PlotPolymerTube(pol,'showSpheres',false,'tubeRadius',.5,'method','spline','colormap',cmap2,'center',false,'lightOn',false); hold on;
PlotSpheres(pol(tads,:),'r',1,'color','r','lightingOn',true); hold on;
set(gca,'color','w'); axis off;


%%
%%
if 0  % prevent accidential overwrite execution
    f = gcf;
    name = ['U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Images\ExampleNonLE_TADs_vEstick.pdf'];
    exportgraphics(f,name,'ContentType','vector');
end
