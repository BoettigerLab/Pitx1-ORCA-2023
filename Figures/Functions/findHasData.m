function hasData = findHasData(Lmaps)
    pitx = 22;
    ra3 = 35;
    ra4 = 46;
    pen = 55;
    nanCrit = size(Lmaps, 1) - 1; % criteria for not having data. Sometimes the ones without data still have a value at diagonal.
    nanMaps = squeeze(sum(isnan(Lmaps), 1));
    hasP = squeeze(nanMaps(pitx, :) < nanCrit)';
    hasE = squeeze(nanMaps(pen, :) < nanCrit)';
    hasB1B2 = squeeze(Lmaps(ra3, ra4, :) < inf);
    hasData = hasP & hasE & hasB1B2;
end
