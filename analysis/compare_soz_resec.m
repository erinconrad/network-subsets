function compare_soz_resec(soz)



contig_text = 'contiguous';
sec_text = 'sec_neg5';


metrics = {'rho_cc','rho_bc','rho_ns','eff','sync'};
dists = {'dist_soz','dist_resec','overlap_soz','overlap_resec'};
    

for dist = 1:length(dists)
    figure
    set(gcf,'Position',[8 138 1239 667]);
    [ha, pos] = tight_subplot(2, 2, [0.07 0.1], [0.1 0.06],[0.11 0.02]);
    count = 0;
    for metric = 1:length(metrics)
        count = count + 1;
        axes(ha(count))
        
        z = [];
        measure_all = [];
        dist_measure_all = [];
        for i = 1:length(soz)
            base = soz(i).(contig_text).(sec_text);
            measure = base.(metrics{metric})';
            dist_measure = base.(dists{dist})';
            rho =  corr(measure,dist_measure);
            measure_all = [measure_all;measure];
            dist_measure_all = [dist_measure_all;dist_measure];
            
            z = [z;atanh(rho)];
        end
        
        comb_z = nanmean(z);
        comb_rho = tanh(comb_z);
        comb_p = normcdf(comb_z);
        
        fprintf('For %s and %s, %d had positive correlation and %d had negative.\n\n',...
            metrics{metric},dists{dist},sum(z>0),sum(z<0));
        
        scatter(dist_measure_all,measure_all,50,'DisplayName',sprintf('%s %s:\nr = %1.2f, p = %1.2e',...
            dists{dist},metrics{metric},comb_rho,comb_p));
        l = legend('location','southeast');
        set(l,'Interpreter', 'none');
        title(sprintf('%s %s',dists{dist},metrics{metric}),'Interpreter', 'none');
        set(gca,'fontsize',20)
        
    end
end
    
   


end