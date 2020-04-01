import yaml

gene_dict = {}

filename = input('Please input filename of desired output: ') + '.yml'
gene_dict.update({'length_of_genome': int(input('Please enter the desired length of the genome: '))})
num_genes = int(input('Please enter desired number of genes: '))
gene_dict.update({'num_genes': num_genes})
for i in range(1, num_genes+1):
    gene_dict.update({'gene{}'.format(i): dict(start=int(input('Please input gene{} starting point: '.format(i))), stop=int(input('Please input gene{} ending point: '.format(i))))})

with open('../../data/gene_parameters/'+filename, 'w') as gene_info:
    yaml.dump(gene_dict, gene_info, default_flow_style=False)

print('Success!')
