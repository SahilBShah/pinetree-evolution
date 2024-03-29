import sys
import yaml

def create_yaml(starting_file, gene_file):
    """
    Creates a YAML file with the same starting point each time that can be edited as simulations occur.
    Input(s):
    starting_file is the file name of the initial taml file containing the termplate genome's architectural information.
    gene_file is the user-inputted file containing the desired gene locations, number of genes, and initial genome length.
    Output(s):
    Outputs the intial yaml file containing the genome's information to the output directory.
    """

    data = dict(
                promoter_0 = dict(start=1, stop=35, current_strength=10e7, previous_strength=0),
                region_0 = dict(start=35, stop=gene_file['gene_1']['start']-30),
    )

    if data['region_0']['start'] > data['region_0']['stop']:
        print('Intragenic region is not long enough before gene 1.')
        sys.exit()

    #Adds promoters to yaml file
    for num_prom in range(1, gene_file['num_genes']):
        if 'promoter_{}'.format(num_prom) not in gene_file.keys():
            data['promoter_{}'.format(num_prom)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
        else:
            data['promoter_{}'.format(num_prom)] = gene_file['promoter_{}'.format(num_prom)]
    #Adds terminators to yaml file
    for num_term in range(1, gene_file['num_genes']+1):
        if 'terminator_{}'.format(num_term) not in gene_file.keys():
            data['terminator_{}'.format(num_term)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
        else:
            data['terminator_{}'.format(num_term)] = gene_file['terminator_{}'.format(num_term)]
    #Adds rnases to yaml file
    for num_rnase in range(gene_file['num_genes']):
        if 'rnase_{}'.format(num_rnase) not in gene_file.keys():
            data['rnase_{}'.format(num_rnase)] = dict(start=0, stop=0, current_strength=0, previous_strength=0)
        else:
            data['rnase_{}'.format(num_rnase)] = gene_file['rnase_{}'.format(num_rnase)]
    #Adds regions to yaml file
    for num_region in range(1, gene_file['num_genes']):
        if 'region_{}'.format(num_region) not in gene_file.keys():
            data['region_{}'.format(num_region)] = dict(start=gene_file['gene_{}'.format(num_region)]['stop'], stop=gene_file['gene_{}'.format(num_region+1)]['start']-30)
        else:
            data['region_{}'.format(num_region)] = gene_file['region_{}'.format(num_region)]
        if num_region == gene_file['num_genes'] - 1:
            data['region_{}'.format(num_region+1)] = dict(start=gene_file['gene_{}'.format(num_region+1)]['stop'], stop=gene_file['length_of_genome'])
        if data['region_{}'.format(num_region)]['start'] > data['region_{}'.format(num_region)]['stop']:
            print('Intragenic region is not long enough between genes ' + str(num_region) + ' and ' + str(num_region+1) + '.')
            sys.exit()
    #Adds genes to yaml file
    for gene_num in range(1, gene_file['num_genes']+1):
        if 'gene_{}'.format(gene_num) not in gene_file.keys():
            data['gene_{}'.format(gene_num)] = dict(start=gene_file['gene_{}'.format(gene_num)]['start'], stop=gene_file['gene_{}'.format(gene_num)]['stop'])
        else:
            data['gene_{}'.format(gene_num)] = gene_file['gene_{}'.format(gene_num)]
    data['length_of_genome'] = gene_file['length_of_genome']
    data['num_genes'] = gene_file['num_genes']
    #Exports yaml file with data
    with open(starting_file, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

    return
