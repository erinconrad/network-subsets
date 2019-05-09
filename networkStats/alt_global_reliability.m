function reliability = alt_global_reliability(std_error,std_true)

variance_error = std_error.^2;
variance_true = repmat(std_true.^2,size(std_error,1),1,size(std_error,3));

reliability = variance_true./(variance_true+variance_error);

end


