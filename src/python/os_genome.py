import os

class call_pt:

    def pt_call(output_dir):

        os.system('./genome_simulator.py {}'.format(output_dir))
