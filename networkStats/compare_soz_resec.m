function compare_soz_resec(soz)

all_sync = [];
all_rho_cc = [];
all_dist_soz = [];
all_eff = [];

contig_text = 'contiguous';
sec_text = 'sec_neg5';

for i = 1:length(soz)
    all_sync = [all_sync,soz(i).(contig_text).(sec_text).sync];
    all_eff = [all_eff,soz(i).(contig_text).(sec_text).eff];
    
    all_rho_cc = [all_rho_cc,soz(i).(contig_text).(sec_text).rho_cc];
    all_dist_soz = [all_dist_soz,soz(i).(contig_text).(sec_text).dist_soz];
    
end

figure
scatter(all_dist_soz,all_rho_cc);

[rho,p] = corr(all_dist_soz(~isnan(all_rho_cc)&~isnan(all_dist_soz))',...
    all_rho_cc(~isnan(all_rho_cc)&~isnan(all_dist_soz))','Type','Spearman')
xlabel('Distance from SOZ');
ylabel('Correlation'); 


figure
scatter(all_dist_soz,all_sync);
[rho,p] = corr(all_dist_soz(~isnan(all_eff)&~isnan(all_dist_soz))',...
    all_eff(~isnan(all_eff)&~isnan(all_dist_soz))','Type','Spearman')
xlabel('Distance from SOZ');
ylabel('Value'); 


end