function [rho,SMC] = doStats(orig,current)

    if sum(isnan(current)) == length(current)
        rho = nan;
        SMC = nan;
    else

        % Spearman rank coefficient
        rho = corr(orig(~isnan(current)),current(~isnan(current)),...
            'Type','Spearman');

        % Get binary versions for Simple matching coefficient
        orig_bin = orig(~isnan(current)) > 0;
        current_bin = current(~isnan(current)) > 0;

        % Get simple matching coefficient
        SMC = simple_matching_coefficient(orig_bin,current_bin);
        
    end
    
end