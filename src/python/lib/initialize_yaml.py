import sys
import yaml

def create_yaml(starting_file, gene_file):
    """
    Creates a YAML file with the same starting point each time that can be edited as simulations occur.
    Input(s):
    starting_file is the file name of the initial taml file containing the termplate genome's architectural information.
    gene_file is the user-inputted file containing the desired gene locations, number of genes, and initial genome length.
    """

    data = dict(
                promoter_0 = dict(start=1, stop=35, current_strength=10e8, previous_strength=0),
                region_0 = dict(start=35, stop=gene_file['gene_1']['start']-30),
    )

    if data['region_0']['start'] > data['region_0']['stop']:
        print('Intergenic region is not long before gene 1.')
        sys.exit()

    #Adds promoters to yaml file
    for num_prom in range(1, gene_file['num_genes']):
        data['promoter_{}'.format(num_prom)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
    #Adds terminators to yaml file
    for num_term in range(1, gene_file['num_genes']+1):
        data['terminator_{}'.format(num_term)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
    #Adds rnases to yaml file
    for num_rnase in range(gene_file['num_genes']):
        data['rnase_{}'.format(num_rnase)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
    #Adds regions to yaml file
    for num_region in range(1, gene_file['num_genes']):
        data['region_{}'.format(num_region)] = dict(start=gene_file['gene_{}'.format(num_region)]['stop'], stop=gene_file['gene_{}'.format(num_region+1)]['start']-30)
        if num_region == gene_file['num_genes'] - 1:
            data['region_{}'.format(num_region+1)] = dict(start=gene_file['gene_{}'.format(num_region+1)]['stop'], stop=gene_file['length_of_genome'])
        if data['region_{}'.format(num_region)]['start'] > data['region_{}'.format(num_region)]['stop']:
            print('Intergenic region is not long enough between genes ' + str(num_region) + ' and ' + str(num_region+1) + '.')
            sys.exit()
    #Adds genes to yaml file
    for gene_num in range(1, gene_file['num_genes']+1):
        data['gene_{}'.format(gene_num)] = dict(start=gene_file['gene_{}'.format(gene_num)]['start'], stop=gene_file['gene_{}'.format(gene_num)]['stop'])
    data['length_of_genome'] = gene_file['length_of_genome']
    data['num_genes'] = gene_file['num_genes']
    #Exports yaml file with data
    with open(starting_file, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

    return
