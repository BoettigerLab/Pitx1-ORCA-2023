function [] = plotSimPolymer(currPoly, ctcf3D, lefs, style, lefMethod)
    figure();
    %Define colormap
    if size(currPoly,1) == 75
        pel = 14;   pitx = 22;  ra3 = 35;   ra4 = 46;   pen = 55;   neurog = 64;
        cmap = zeros(75,3) + [120 120 120]/255;
        cmap(pitx:ra3,:) = repmat([255 127 39]/255,ra3-pitx+1,1);
        cmap(ra3+1:ra4,:) = repmat([255 230 80]/255,ra4-ra3,1);
        cmap(ra4+1:58,:) = repmat([89 216 52]/255,58-ra4,1);
    else
        segLen = floor(size(currPoly,1)/5);
        pitx = segLen;
        ra3 = segLen*2;
        ra4 = segLen*3;
        pen = segLen*4;
        cmap = zeros(size(currPoly,1),3) + [120 120 120]/255;
        cmap(pitx+1:ra3,:) = repmat([255 127 39]/255,ra3-pitx,1);
        cmap(ra3+1:ra4,:) = repmat([255 230 80]/255,ra4-ra3,1);
        cmap(ra4+1:pen,:) = repmat([89 216 52]/255,pen-ra4,1);
    end
    
    if strcmp(style, "thick")
        tubeRadius = 0.5;
        sphereRadius = 1;
    elseif strcmp(style, "thin")
        tubeRadius = 0.15;
        sphereRadius = 0.3;
    end
    PlotPolymerTube_TzuChiao(currPoly,'tubeRadius',tubeRadius,'showSpheres',false,...
        'colormap',cmap,'alpha',.7,'method','pchip','center',false,'lightOn',false);
    set(gca,'color','w'); hold on;
    if size(currPoly,1) == 75
        PlotSpheres(ctcf3D,'color',[91 199 206]/255,'r',sphereRadius, 'lightingOn', false); hold on;
    else
        PlotSpheres(ctcf3D,'color',[255 127 39; 255 230 80; 180 210 40; 89 216 52]/255,'r',sphereRadius, 'lightingOn', false); hold on;
    end
    PlotLEFs(lefs,'color',[0 0 0],'r',sphereRadius,'alpha',.5, 'method', lefMethod); hold on;
    axis equal;
end