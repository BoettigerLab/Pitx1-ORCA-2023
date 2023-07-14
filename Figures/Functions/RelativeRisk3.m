function outCell = RelativeRisk3(condition, exposure, has)
    % condition: outcome
    % exposure: predictor
    % has: indices that should be taken into account
    % outCell: [rr, CI lower bound, CI upper bound, N]
    outCell = zeros(4,1);

    A = sum(condition & exposure & has);
    B = sum(~condition & exposure & has);
    C = sum(condition & ~exposure & has);
    D = sum(~condition & ~exposure & has);
    rr = (A/(A+B))/(C/(C+D));
    
    cI = 0.95;
    Z = norminv([(1 - cI)/2, 0.5 + 0.5*cI]);
    Z = Z(2);

    logCI = Z*((B/A)/(A+B) + (D/C)/(C+D))^0.5;
    rr_CI = [exp(log(rr) - logCI), exp(log(rr) + logCI)];
    outCell(1) = rr;
    outCell(2:3) = rr_CI;
    outCell(4) = sum(has);
end
