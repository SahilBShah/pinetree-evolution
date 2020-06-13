import pinetree as pt
import random
import yaml


def pt_call(output_dir, genome_tracker_new):

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(random.randint(0, 10e6))
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=35)
    sim.add_ribosome(copy_number=100, speed=30, footprint=30)

    plasmid = pt.Genome(name="plasmid", length=genome_tracker_new['length_of_genome'],
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)
    #Promoters
    for i in range(genome_tracker_new['num_genes']):
        if genome_tracker_new['promoter_{}'.format(i)]['start'] > 0:
            plasmid.add_promoter(name="p{}".format(i), start=genome_tracker_new['promoter_{}'.format(i)]['start'], stop=genome_tracker_new['promoter_{}'.format(i)]['stop'],
                                 interactions={"rnapol": genome_tracker_new['promoter_{}'.format(i)]['current_strength']})

    #RNases
    for i in range(genome_tracker_new['num_genes']):
        if genome_tracker_new['rnase_{}'.format(i)]['start'] > 0:
            plasmid.add_rnase_site(start=genome_tracker_new['rnase_{}'.format(i)]['start'], stop=genome_tracker_new['rnase_{}'.format(i)]['stop'])

    #Terminators
    for i in range(1, genome_tracker_new['num_genes']+1):
        if genome_tracker_new['terminator_{}'.format(i)]['start'] > 0:
            plasmid.add_terminator(name="t{}".format(i), start=genome_tracker_new['terminator_{}'.format(i)]['start'], stop=genome_tracker_new['terminator_{}'.format(i)]['stop'],
                                   efficiency={"rnapol": genome_tracker_new['terminator_{}'.format(i)]['current_strength']})

    #Genes
    for i in range(1, genome_tracker_new['num_genes']+1):
        if genome_tracker_new['gene_{}'.format(i)]['start'] > 0:
            plasmid.add_gene(name="protein{}".format(i), start=genome_tracker_new['gene_{}'.format(i)]['start'], stop=genome_tracker_new['gene_{}'.format(i)]['stop'],
                             rbs_start=(genome_tracker_new['gene_{}'.format(i)]['start']-30), rbs_stop=genome_tracker_new['gene_{}'.format(i)]['start']-1, rbs_strength=1e7)

    sim.register_genome(plasmid)
    sim.simulate(time_limit=250, time_step=1, output=output_dir+'expression_pattern.tsv')

def pt_call_alt(output_dir, genome_tracker_new):

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(random.randint(0, 10e6))
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=35)
    sim.add_ribosome(copy_number=100, speed=30, footprint=30)

    plasmid = pt.Genome(name="plasmid", length=genome_tracker_new['length_of_genome'],
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)
    #Promoters
    for i in range(genome_tracker_new['num_genes']):
        if genome_tracker_new['promoter_{}'.format(i)]['start'] > 0:
            plasmid.add_promoter(name="p{}".format(i), start=genome_tracker_new['promoter_{}'.format(i)]['start'], stop=genome_tracker_new['promoter_{}'.format(i)]['stop'],
                                 interactions={"rnapol": genome_tracker_new['promoter_{}'.format(i)]['current_strength']})

    #RNases
    for i in range(genome_tracker_new['num_genes']):
        if genome_tracker_new['rnase_{}'.format(i)]['start'] > 0:
            plasmid.add_rnase_site(name='r{}'.format(i), start=genome_tracker_new['rnase_{}'.format(i)]['start'], stop=genome_tracker_new['rnase_{}'.format(i)]['stop'], rate=genome_tracker_new['rnase_{}'.format(i)]['current_strength'])

    #Terminators
    for i in range(1, genome_tracker_new['num_genes']+1):
        if genome_tracker_new['terminator_{}'.format(i)]['start'] > 0:
            plasmid.add_terminator(name="t{}".format(i), start=genome_tracker_new['terminator_{}'.format(i)]['start'], stop=genome_tracker_new['terminator_{}'.format(i)]['stop'],
                                   efficiency={"rnapol": genome_tracker_new['terminator_{}'.format(i)]['current_strength']})

    #Genes
    for i in range(1, genome_tracker_new['num_genes']+1):
        if genome_tracker_new['gene_{}'.format(i)]['start'] > 0:
            plasmid.add_gene(name="protein{}".format(i), start=genome_tracker_new['gene_{}'.format(i)]['start'], stop=genome_tracker_new['gene_{}'.format(i)]['stop'],
                             rbs_start=(genome_tracker_new['gene_{}'.format(i)]['start']-30), rbs_stop=genome_tracker_new['gene_{}'.format(i)]['start']-1, rbs_strength=1e7)

    sim.register_genome(plasmid)
    sim.simulate(time_limit=250, time_step=1, output=output_dir+'expression_pattern.tsv')

    return
