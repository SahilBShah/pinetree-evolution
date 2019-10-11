#!/usr/bin/env python
import argparse
import pinetree as pt
import random
import yaml


parser = argparse.ArgumentParser()
parser.add_argument('output_directory')
args = parser.parse_args()
output_dir = args.output_directory

with open(output_dir+'new_gene.yml') as f:
    genome_tracker_new = yaml.safe_load(f)

sim = pt.Model(cell_volume=8e-16)
sim.seed(random.randint(0, 10e6))
sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
sim.add_ribosome(copy_number=100, speed=30, footprint=10)

plasmid = pt.Genome(name="plasmid", length=450,
                    transcript_degradation_rate=1e-2,
                    transcript_degradation_rate_ext=1e-2,
                    rnase_speed=20,
                    rnase_footprint=10)
plasmid.add_promoter(name="p1", start=1, stop=10,
                     interactions={"rnapol": genome_tracker_new['promoter1']['current_strength']})
if genome_tracker_new['promoter2']['start'] > 0:
    plasmid.add_promoter(name="p2", start=genome_tracker_new['promoter2']['start'], stop=genome_tracker_new['promoter2']['stop'],
                         interactions={"rnapol": genome_tracker_new['promoter2']['current_strength']})
if genome_tracker_new['promoter3']['start'] > 0:
    plasmid.add_promoter(name="p3", start=genome_tracker_new['promoter3']['start'], stop=genome_tracker_new['promoter3']['stop'],
                         interactions={"rnapol": genome_tracker_new['promoter3']['current_strength']})
if genome_tracker_new['rnase1']['start'] > 0:
    plasmid.add_rnase_site(start=genome_tracker_new['rnase1']['start'], stop=genome_tracker_new['rnase1']['stop'])
if genome_tracker_new['rnase2']['start'] > 0:
    plasmid.add_rnase_site(start=genome_tracker_new['rnase2']['start'], stop=genome_tracker_new['rnase2']['stop'])
if genome_tracker_new['rnase3']['start'] > 0:
    plasmid.add_rnase_site(start=genome_tracker_new['rnase3']['start'], stop=genome_tracker_new['rnase3']['stop'])
if genome_tracker_new['terminator1']['start'] > 0:
    plasmid.add_terminator(name="t1", start=genome_tracker_new['terminator1']['start'], stop=genome_tracker_new['terminator1']['stop'],
                           efficiency={"rnapol": genome_tracker_new['terminator1']['current_strength']})
if genome_tracker_new['terminator2']['start'] > 0:
    plasmid.add_terminator(name="t2", start=genome_tracker_new['terminator2']['start'], stop=genome_tracker_new['terminator2']['stop'],
                           efficiency={"rnapol": genome_tracker_new['terminator2']['current_strength']})
if genome_tracker_new['terminator3']['start'] > 0:
    plasmid.add_terminator(name="t3", start=genome_tracker_new['terminator3']['start'], stop=genome_tracker_new['terminator3']['stop'],
                           efficiency={"rnapol": genome_tracker_new['terminator3']['current_strength']})
plasmid.add_gene(name="proteinX", start=genome_tracker_new['geneX']['start'], stop=genome_tracker_new['geneX']['stop'],
                 rbs_start=(genome_tracker_new['geneX']['start']-15), rbs_stop=genome_tracker_new['geneX']['start'], rbs_strength=1e7)
plasmid.add_gene(name="proteinY", start=genome_tracker_new['geneY']['start'], stop=genome_tracker_new['geneY']['stop'],
                 rbs_start=(genome_tracker_new['geneY']['start']-15), rbs_stop=genome_tracker_new['geneY']['start'], rbs_strength=1e7)
plasmid.add_gene(name="proteinZ", start=genome_tracker_new['geneZ']['start'], stop=genome_tracker_new['geneZ']['stop'],
                 rbs_start=(genome_tracker_new['geneZ']['start']-15), rbs_stop=genome_tracker_new['geneZ']['start'], rbs_strength=1e7)
sim.register_genome(plasmid)
sim.simulate(time_limit=240, time_step=1, output=output_dir+'three_genes_replicated.tsv')
