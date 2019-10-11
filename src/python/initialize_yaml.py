import yaml

def create_yaml(starting_file):
    """
    Creates a YAML file with the same starting point each time that can be edited as simulations occur
    """

    data = dict(
                geneX = dict(start=26, stop=121),
                geneY = dict(start=159, stop=280),
                geneZ = dict(start=319, stop=449),
                promoter1 = dict(start=1, stop=10, current_strength=10e8, previous_strength=0),
                promoter2 = dict(start=0, stop=0, current_strength=0, previous_strength=0),
                promoter3 = dict(start=0, stop=0, current_strength=0, previous_strength=0),
                rnase1 = dict(start=0, stop=0),
                rnase2 = dict(start=0, stop=0),
                rnase3 = dict(start=0, stop=0),
                terminator1 = dict(start=0, stop=0, current_strength=0, previous_strength=0),
                terminator2 = dict(start=0, stop=0, current_strength=0, previous_strength=0),
                terminator3 = dict(start=0, stop=0, current_strength=0, previous_strength=0),
                region1 = dict(start=11, stop=25),
                region2a = dict(start=122, stop=132),
                region2b = dict(start=133, stop=143),
                region2c = dict(start=144, stop=159),
                region3a = dict(start=281, stop=291),
                region3b = dict(start=292, stop=302),
                region3c = dict(start=303, stop=318),
    )

    with open(starting_file, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
    outfile.close()
