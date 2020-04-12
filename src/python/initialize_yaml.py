import yaml

def create_yaml(starting_file, gene_file):
    """
    Creates a YAML file with the same starting point each time that can be edited as simulations occur.
    """

    data = dict(
                promoter_0 = dict(start=1, stop=10, current_strength=10e8, previous_strength=0),
                region_0 = dict(start=10, stop=gene_file['gene_1']['start']),
    )
    #Adds promoters to yaml file
    for num_prom in range(1, gene_file['num_genes']):
        data.update({'promoter_{}'.format(num_prom): dict(start=0, stop=0, current_strength=0, previous_strength=0)})
    #Adds terminators to yaml file
    for num_term in range(1, gene_file['num_genes']+1):
        data.update({'terminator_{}'.format(num_term): dict(start=0, stop=0, current_strength=0, previous_strength=0)})
    #Adds rnases to yaml file
    for num_rnase in range(gene_file['num_genes']):
        data.update({'rnase_{}'.format(num_rnase): dict(start=0, stop=0, current_strength=0, previous_strength=0)})
    #Adds regions to yaml file
    for num_region in range(1, gene_file['num_genes']):
        data.update({'region_{}'.format(num_region): dict(start=gene_file['gene_{}'.format(num_region)]['stop'], stop=gene_file['gene_{}'.format(num_region+1)]['start'])})
        if num_region == gene_file['num_genes'] - 1:
            data.update({'region_{}'.format(num_region+1): dict(start=gene_file['gene_{}'.format(num_region+1)]['stop'], stop=gene_file['length_of_genome'])})
    #Adds genes to yaml file
    for gene_num in range(1, gene_file['num_genes']+1):
        data.update({'gene_{}'.format(gene_num): dict(start=gene_file['gene_{}'.format(gene_num)]['start'], stop=gene_file['gene_{}'.format(gene_num)]['stop'])})
    data.update({'length_of_genome': gene_file['length_of_genome']})
    data.update({'num_genes': gene_file['num_genes']})
    #Exports yaml file with data
    with open(starting_file, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

    return
