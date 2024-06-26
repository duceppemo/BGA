import os
import sys
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from bga_methods import Methods
from nanopore_methods import NanoporeMethods
from illumina_methods import IlluminaMethods
from assembly_qc import AssemblyQcMethods
from importlib import resources


__author__ = 'duceppemo'
__version__ = '0.1'


"""
# Create virtual environment
mamba create -n BGA -c conda-forge -c bioconda -c plotly -y python=3.10.8 nextpolish=1.4.1 bwa=0.7.17 samtools=1.17 \
    porechop=0.2.4 filtlong=0.2.1 minimap2=2.26 flye=2.9.2 shasta=0.11.1 qualimap=2.2.2d bbmap=39.01 bandage=0.8.1 \
    fastp=0.22.0 ntedit=1.3.5 polypolish=0.5.0 pandas=1.5.3 seqtk=1.4 quast=5.2.0 medaka=1.8.0 mummer4=4.0.0rc1 \
    gnuplot=5.4.5 plotly=5.15.0 pigz=2.6 falco=1.2.1

plotly=5.14.1
r-base=4.2.2 r-optparse=1.7.3 r-ggplot2=3.4.2 r-plotly=4.10.1

pip install mummer-idotplot

# Activate virtual environment
conda activate BGA

"""


class BGA(object):
    def __init__(self, args):
        # I/O
        self.long_reads = os.path.abspath(args.long_reads)
        self.short_reads = ''
        if args.short_reads:
            self.short_reads = os.path.abspath(args.short_reads)
        self.output_folder = os.path.abspath(args.output)
        self.ref_size = args.size
        self.assembler = args.assembler
        self.min_size = args.min_size
        self.model = args.model

        # Flags
        self.trim_long = args.trim_long
        self.filter_long = args.filter_long
        self.trim_short = args.trim_short
        self.polish = args.polish

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Resources
        # self.rscript = resources.Resource(os.path('dependencies', 'mummerCoordsDotPlotly.R'))
        # Data
        self.sample_dict = dict()

        # Log
        self.log_file = self.output_folder + '/log.txt'

        # Run
        os.chdir(self.output_folder)
        self.run()

    def run(self):
        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        # Check input file compatibility
        Methods.check_input(self.long_reads)
        if self.short_reads:
            Methods.check_input(self.short_reads)

        Methods.check_dependencies(self.output_folder)
        print('\tAll good!')

        ############################################################

        # Step completion report files
        done_trimming_long = self.output_folder + '/done_trimming_long'
        done_filtering = self.output_folder + '/done_filtering_long'
        done_assembling = self.output_folder + '/done_assembling_long'
        done_polishing_long = self.output_folder + '/done_polishing_long'
        done_trimming_short = self.output_folder + '/done_trimming_short'
        done_polishing_short = self.output_folder + '/done_polishing_short'
        done_qc = self.output_folder + '/done_qc'

        # Output folders to create
        trimmed_folder = self.output_folder + '/1_trimmed_long/'
        filtered_folder = self.output_folder + '/2_filtered_long/'
        assembled_folder = self.output_folder + '/3_assembled_long/'
        polished_long_folder = self.output_folder + '/4_polished_long/'
        illumina_trimmed_folder = self.output_folder + '/5a_trimmed_short/'
        illumina_trimmed_report_folder = self.output_folder + '/5b_trimmed_short_report/'
        polished_short_folder = self.output_folder + '/6_polished_short/'
        assembly_qc_folder = self.output_folder + '/7_assembly_qc/'
        quast_folder = assembly_qc_folder + 'quast/'
        qualimap_long_folder = assembly_qc_folder + 'qualimap_long_reads/'
        qualimap_short_folder = assembly_qc_folder + 'qualimap_short_reads/'
        mummer_long_folder = assembly_qc_folder + '/mummer_long/'
        mummer_short_folder = assembly_qc_folder + '/mummer_short/'
        gfa_folder = assembly_qc_folder + 'assembly_graphs/'
        coverage_long_folder = assembly_qc_folder + '/coverage_long/'
        coverage_short_folder = assembly_qc_folder + '/coverage_short/'

        # Get input files and place info in dictionary
        if self.long_reads:
            # Get fastq files
            self.sample_dict = Methods.get_nanopore_files(self.long_reads)

            # Drop the unclassified if present
            # We don't want to assemble those reads because we don't know from which barcode they come from
            self.sample_dict = {k: v for (k, v) in self.sample_dict.items() if 'unclassified' not in k}

            # Trim
            if self.trim_long:
                print('Trimming long reads with Porechop...')
                NanoporeMethods.run_porechop_parallel(self.sample_dict, trimmed_folder,
                                                      self.cpu, self.parallel, done_trimming_long)

            # Filter
            if self.filter_long:
                print('Filtering long reads with Filtlong...')
                NanoporeMethods.run_filtlong_parallel(self.sample_dict, filtered_folder,
                                                      self.ref_size, self.parallel, done_filtering)

            # Assemble
            if self.assembler == 'flye':
                print('Assembling long reads with Flye...')
                NanoporeMethods.assemble_flye_parallel(self.sample_dict, assembled_folder, gfa_folder, self.ref_size,
                                                       self.min_size, self.cpu, self.parallel, done_assembling,
                                                       self.output_folder)
            else:  # elif self.assembler == 'shasta':
                print('Assembling long reads with Shasta...')
                # sample_dict, output_folder, gfa_folder, min_size, cpu, parallel, flag
                NanoporeMethods.assemble_shasta_parallel(self.sample_dict, assembled_folder, gfa_folder,
                                                         self.min_size, self.cpu, self.parallel, done_assembling,
                                                         self.output_folder)

            # Long read polishing
            print('Long read polishing with Medaka...')
            # NanoporeMethods.polish_medaka_parallel(self.sample_dict, polished_long_folder, self.model,
            #                                        self.cpu, self.parallel, done_polishing_long)
            NanoporeMethods.polish_medaka_loop(self.sample_dict, polished_long_folder, self.model,
                                               self.cpu, done_polishing_long)
        else:
            raise Exception('You must provide long reads for the assembly.')

        if self.polish:
            if self.short_reads:
                # Get fastq files
                self.sample_dict = Methods.get_illumina_files(self.short_reads, self.sample_dict)

                # Trim
                if self.trim_short:
                    print('Trimming short reads with FastP...')
                    IlluminaMethods.trim_illumina_fastp_paired_parallel(self.sample_dict, illumina_trimmed_folder,
                                                                        illumina_trimmed_report_folder,
                                                                        self.cpu, self.parallel, done_trimming_short)

                # Polish
                print('Short read polishing with Polypolish and pyPOLCA...')
                IlluminaMethods.polish_parallel(self.sample_dict, polished_short_folder,
                                                self.cpu, self.parallel, done_polishing_short)
            else:
                raise Exception('You must provide paired-end short read data in order to perform short read polishing.')

        # QC
        print('Performing assembly QC...')

        # Pre- / post-long read polishing comparison
        # sample_dict, ref_folder, query_folder, output_folder, cpu, parallel
        print('\tRunning Nucmer - raw VS long read polished')
        # AssemblyQcMethods.run_last_parallel(self.sample_dict, assembled_folder, polished_long_folder,
        #                                     assembly_qc_folder, self.cpu, self.parallel)
        AssemblyQcMethods.run_nucmer_medaka_parallel(self.sample_dict, mummer_long_folder, self.cpu, self.parallel)

        # Qualimap long reads
        print('\tRunning Qualimap - long reads')
        AssemblyQcMethods.map_minimap2_parallel(self.sample_dict, qualimap_long_folder, 'nanopore',
                                                self.cpu, self.parallel)

        AssemblyQcMethods.run_qualimap_parallel(self.sample_dict, qualimap_long_folder,
                                                self.cpu, self.mem, self.parallel)

        # QUAST
        print('\tRunning QUAST')
        AssemblyQcMethods.run_quast(self.sample_dict, quast_folder, self.cpu)

        if self.polish:
            print('\tRunning Nucmer - long read polished VS short read polished')
            # AssemblyQcMethods.run_last_parallel(self.sample_dict, assembled_folder, polished_long_folder,
            #                                     assembly_qc_folder, self.cpu, self.parallel)
            AssemblyQcMethods.run_nucmer_polypolish_parallel(self.sample_dict, mummer_short_folder,
                                                             self.cpu, self.parallel)

            # Qualimap short reads
            print('\tRunning Qualimap - short reads')
            AssemblyQcMethods.map_minimap2_parallel(self.sample_dict, qualimap_long_folder, 'illumina',
                                                    self.cpu, self.parallel)

            AssemblyQcMethods.run_qualimap_parallel(self.sample_dict, qualimap_short_folder, self.cpu, self.mem,
                                                    self.parallel)

            print('\tComputing short read coverage')
            AssemblyQcMethods.short_read_coverage_parallel(self.sample_dict, coverage_short_folder, self.cpu,
                                                           self.parallel)

        print('\tComputing long read coverage')
        AssemblyQcMethods.long_read_coverage_parallel(self.sample_dict, coverage_long_folder, self.cpu,
                                                      self.parallel)

        print('Done!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Long read or hybrid bacterial genome assembly pipeline.')
    parser.add_argument('-l', '--long-reads', metavar='/path/to/nanopore_folder',
                        required=True, type=str,
                        help='Folder that contains the fastq files from Nanopore. Mandatory.')
    parser.add_argument('-s', '--short-reads', metavar='/path/to/illumina_folder',
                        required=False, type=str,
                        help='Folder that contains the fastq files from Illumina. Illumina file names must match '
                             'Nanopore files names (everything before the first underscore ("_"). Optional.')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-a', '--assembler',
                        required=False, default='flye',
                        choices=['flye', 'shasta'],
                        type=str,
                        help='Assembly method. Default "flye". Optional.')
    parser.add_argument('--min-size', metavar='3000',
                        required=False,
                        type=int,
                        help='Minimum read size for Shasta assembler or minimum read overlap for Flye. '
                             'Shasta default is 3,000. Flye default is "auto". Optional.')
    parser.add_argument('--size', metavar='5000000',
                        required=False,
                        type=int,
                        help='Override automatically detected reference size for Flye. If entered manually, '
                             'it will also be used to during the long read filtering step (to retain the 100x to reads,'
                             ' if coverage allows. Optional.')
    parser.add_argument('--trim-long',
                        required=False, action='store_true',
                        help='Trim long reads with Porechop prior assembly. Default is False.')
    parser.add_argument('--filter-long',
                        required=False, action='store_true',
                        help='Filter long reads with Filtlong prior assembly. Drop bottom 5%%. Default is False.')
    parser.add_argument('--model',
                        type=str, required=False, default='r1041_e82_400bps_sup_g615',
                        help='Medaka model. Default is for R10.4.1 v14 flowcell.')
    parser.add_argument('--show-models',
                        action='help',
                        help='Show available models for Medaka polishing and exit.')
    parser.add_argument('--polish',
                        required=False, action='store_true',
                        help='Polish long read assembly with Illumina paired-end data. Default is False.')
    parser.add_argument('--trim-short',
                        required=False, action='store_true',
                        help='Trim paired-end Illumina reads with FastP prior polishing. Default is False.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Keep low if your computer has low memory. '
                             'Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    if len(sys.argv) == 2 and sys.argv[1] == '--show-models':
        medaka_model_list = NanoporeMethods.print_medaka_models()
        print('Available medaka models:\n{}'.format(', '.join(medaka_model_list)))
        sys.exit()

    # Get the arguments into an object
    arguments = parser.parse_args()

    BGA(arguments)
