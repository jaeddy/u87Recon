function biomass = read_biomass(model, biomass_file)

    %% read biomass file
    fid = fopen(biomass_file, 'r');
    t = textscan(fid, '%f %s');
    fclose(fid);
    coefficients = t{1};
    metabolites = t{2};
    metabolites = regexprep(metabolites,'\-','_');
    metabolites = regexprep(metabolites,'DASH_','');
    
    % fill in biomass vector
    biomass = zeros(size(model.mets));
    for i=1:length(metabolites)
        metabolite = char(metabolites(i));
        coefficient = coefficients(i);
        if (coefficient == 0)
            continue;
        end
        metabolite_index = strcmpi(metabolite, model.mets);
%         fprintf('%s\n',metabolite);
%         p(i) = find(metabolite_index);
        metabolite_index_count = nnz(metabolite_index);
        if metabolite_index_count ~= 1
            error('Metabolite %s found %d times', metabolite, metabolite_index_count);
        end
        biomass(metabolite_index) = coefficient;
    end
end
