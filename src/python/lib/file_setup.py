import pandas as pd

#Removes unnecessary rows and columns in produced file
def rearrange_file(file, max_time, num_genes):
    """
    Rearranges files so that the target input file can be compared to the simulated output file so that sum of squares can be easily calculated.
    Input(s):
    file is the transcript abundances file outputted from pinetree.
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    Output(s):
    file is a dataframe with a modified layout.
    """

    protein_species = []

    #File is rearranged to have columns represent each gene's transcript abundances and each row represent the time
    for gene_number in range(1, num_genes+1):
        protein_species.append('protein{}'.format(gene_number))
    file = file[file['species'].isin(protein_species)]
    file = file[['time', 'species', 'transcript']]
    #Time is rounded so that files can be compared to one another
    file['time'] = file['time'].round().astype(int)
    file = file.pivot(index='time', columns='species', values='transcript')
    file = file.fillna(0.0)
    #If a protein is not present at a certain time period then no (0) transcripts are produced
    for protein in protein_species:
        if protein not in file.columns:
            file[protein] = 0.0
    #If none of the proteins are present at a certain time point then no (0) transcripts are produced
    for t in range(1, max_time+1):
        if t not in file.index:
            new_row = pd.Series(data={'protein1': 0.0, 'protein2': 0.0, 'protein3': 0.0}, name=t)
            file = file.append(new_row, ignore_index=False)
    file = file.sort_index()
    file = file[protein_species]

    return file
