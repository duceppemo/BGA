import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from bga_methods import Methods
from nanopore_methods import NanoporeMethods
from illumina_methods import IlluminaMethods


__author__ = 'duceppemo'
__version__ = '0.1'


"""
# Create virtual environment
mamba create -n BGA -c bioconda -y python=3.10.11 nextpolish=1.4.1 bwa=0.7.17 samtools=1.17 porechop=0.2.4 \
    filtlong=0.2.1 minimap2=2.26 flye=2.9.2 shasta=0.11.1 quast=5.2.0 qualimap=2.2.2d bbmap=39.01 bandage=0.8.1 \
    fastp=0.23.2 ntedit=1.3.5 polypolish=0.5.0 pandas=1.5.3 seqtk=1.3 medaka=1.8.0

# Activate virtual environment
conda activate BGA

"""


class BGA(object):
    def __init__(self, args):
        # I/O
        self.long_reads = os.path.abspath(args.long_reads)
        self.short_reads = os.path.abspath(args.short_reads)
        self.output_folder = os.path.abspath(args.output)
        self.ref_size = args.size
        self.assembler = args.assembler
        self.min_size = args.min_size

        # Flags
        self.trim_long = args.trim_long
        self.filter_long = args.filter_long
        self. trim_short = args.trim_short
        self.polish = args.polish

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Data
        self.sample_dict = dict()

        # Run
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

        print('\tAll good!')

        ############################################################

        # Step completion report files
        done_trimming_long = self.output_folder + '/done_trimming_long'
        done_filtering = self.output_folder + '/done_filtering_long'
        done_assembling = self.output_folder + '/done_assembling_long'
        done_trimming_short = self.output_folder + '/done_trimming_short'
        done_polishing = self.output_folder + '/done_polishing'

        # Output folders to create
        trimmed_folder = self.output_folder + '/1_trimmed_long/'
        filtered_folder = self.output_folder + '/2_filtered_long/'
        assembled_folder = self.output_folder + '/3_assembled_long/'
        illumina_trimmed_folder = self.output_folder + '/4a_trimmed_short/'
        illumina_trimmed_report_folder = self.output_folder + '/4b_trimmed_short/'
        polished_folder = self.output_folder + '/5_polished_assemblies/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        # Get input files and place info in dictionary
        if self.long_reads:
            # Get fastq files
            self.sample_dict = Methods.get_nanopore_files(self.long_reads)

            # Drop the unclassified if present
            # We don't want to assemble those reads because we don't know from which barcode they come from
            self.sample_dict = {k: v for (k, v) in self.sample_dict.items() if 'unclassified' not in k}

            # Trim
            if self.trim_long:
                print('Trimming Nanopore reads with Porechop...')
                NanoporeMethods.run_porechop_parallel(self.sample_dict, trimmed_folder,
                                                      self.cpu, self.parallel, done_trimming_long)

            # Filter
            if self.filter_long:
                print('Filtering Nanopore reads with Filtlong...')
                NanoporeMethods.run_filtlong_parallel(self.sample_dict, filtered_folder,
                                                      self.ref_size, self.parallel, done_filtering)

            # Assemble
            if self.assembler == 'flye':
                print('Assembling Nanopore reads with Flye...')
                NanoporeMethods.assemble_flye_parallel(self.sample_dict, assembled_folder, self.ref_size,
                                                       self.min_size, self.cpu, self.parallel, done_assembling)
                NanoporeMethods.flye_assembly_stats(assembled_folder, self.output_folder)
            else:  # elif self.assembler == 'shasta':
                print('Assembling Nanopore reads with Shasta...')
                NanoporeMethods.assemble_shasta_parallel(self.sample_dict, assembled_folder,
                                                         self.ref_size, self.cpu, self.parallel, done_assembling)

        else:
            raise Exception('You must provide long reads from Nanopore.')

        if self.short_reads:
            # Get fastq files
            self.sample_dict = Methods.get_illumina_files(self.short_reads, self.sample_dict)

            # Check that we have paired-end reads
            for sample, info_obj in self.sample_dict.items():
                if not info_obj.illumina.raw.r1 and info_obj.illumina.raw.r1:
                    raise Exception('Illumina data must be paired-end (R1 and R2 files required).')

            # Trim
            if self.trim_short:
                print('Trimming short reads with FastP...')
                IlluminaMethods.trim_illumina_fastp_paired_parallel(self.sample_dict, illumina_trimmed_folder,
                                                                    illumina_trimmed_report_folder,
                                                                    self.cpu, self.parallel, done_trimming_short)

            # Polish
            if self.polish:
                print('Polishing long read assembly with short reads...')
                IlluminaMethods.polish(self.sample_dict, polished_folder, self.cpu, done_polishing)
        else:
            raise Exception('You must provide Illumina paired-end data in order to perform polishing.')

        print('\tDone!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Bacterial Genome assembly.')
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
                             'Note that Shasta uses a min read size of 10,000bp by default and Flye automatically '
                             'sets this values based on read length distribution. It is recommended, to set this value '
                             'to something that reflects your read size distribution if using Shasta, say 3,000. '
                             'Optional.')
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
    parser.add_argument('-m', '--model', metavar=)
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

    # Get the arguments into an object
    arguments = parser.parse_args()

    BGA(arguments)
