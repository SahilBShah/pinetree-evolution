import pandas as pd

#Removes unnecessary rows and columns in produced file
def rearrange_file(file):
    """
    Rearranges files so that the target input file can be compared to the simulated output file so that sum of squares can be easily calculated
    """
    protein_species = ['proteinX', 'proteinY', 'proteinZ']
    file = file[file['species'].isin(protein_species)]
    file = file[['time', 'species', 'transcript']]
    file['time'] = file['time'].round().astype(int)
    file = file.pivot(index='time', columns='species', values='transcript')
    file = file.fillna(0.0)
    for protein in protein_species:
        if protein not in file.columns:
            file[protein] = 0.0
    file = file[protein_species]
    return file
