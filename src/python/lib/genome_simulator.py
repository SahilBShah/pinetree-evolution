#common imports
import os
import random
import sys
import time
import yaml

#dependencies
import pinetree as pt


class SupressOutput(object):
    """
    Class used to grab standard output or another stream.
    Adapted from Devan Williams and original code can be found at https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable.
    """

    escape_char = "\b"

    def __init__(self):

        self.origstream = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured
        self.pipe_out, self.pipe_in = os.pipe()

    def __enter__(self):

        self.start()
        return self

    def __exit__(self, type, value, traceback):

        self.stop()

    def start(self):
        """
        Start capturing the stream data.
        """

        self.capturedtext = ""
        #Save a copy of the stream
        self.streamfd = os.dup(self.origstreamfd)
        #Replace the original stream with our write pipe
        os.dup2(self.pipe_in, self.origstreamfd)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """

        #Print the escape character to make the readout method stop
        self.origstream.write(self.escape_char)
        #Flush the stream to make sure all our data goes in before
        self.origstream.flush()
        self.readout()
        #Close the pipe
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        #Restore the original stream
        os.dup2(self.streamfd, self.origstreamfd)
        #Close the duplicate stream
        os.close(self.streamfd)

    def readout(self):
        """
        Read the stream data (one byte at a time) and save the text in `capturedtext`.
        """

        while True:
            char = os.read(self.pipe_out,1).decode(self.origstream.encoding)
            if not char or self.escape_char in char:
                break
            self.capturedtext += char


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

    #Run pinetree simulation
    sim.register_genome(plasmid)
    supress = SupressOutput()
    with supress:
        sim.simulate(time_limit=250, time_step=1, output=output_dir+'expression_pattern.tsv')

    return

def pt_call_alt(output_dir, genome_tracker_new):
    """
    pinetree python interface containing all the information needed to conduct a simulation.
    Input(s):
    output_dir is the path to the directory in which all the saved files are stored by the program.
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    """

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

    #Run pinetree simulation
    sim.register_genome(plasmid)
    supress = SupressOutput()
    with supress:
        sim.simulate(time_limit=250, time_step=1, output=output_dir+'expression_pattern.tsv')

    return
