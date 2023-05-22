import subprocess
import os
from concurrent import futures
from shutil import move
import gzip
from itertools import groupby
from glob import glob
import warnings
import pandas as pd
from bga_methods import Methods


class Sample(object):
    def __getattr__(self, name):
        self.__dict__[name] = Sample()
        return self.__dict__[name]


class NanoporeMethods(object):
    @staticmethod
    def run_porechop(sample, info_obj, trimmed_folder, cpu, flag):
        input_fastq = info_obj.nanopore.raw
        trimmed_fastq = trimmed_folder + '/' + sample + '.fastq.gz'

        if not os.path.exists(flag):
            cmd = ['porechop',
                   '-i', input_fastq,
                   '-o', trimmed_fastq,
                   '--threads', str(cpu),
                   '--check_reads', str(1000)]

            print('\t{}'.format(sample))
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        return sample, trimmed_fastq

    @staticmethod
    def run_porechop_parallel(sample_dict, output_folder, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        # Will skip the actual trimming if flag file exists
        # Still need to run it to update the dictionary with trimmed files
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, int(cpu / parallel), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.run_porechop(*x), args):
                sample_dict[results[0]].nanopore.trimmed = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping filtering long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def run_filtlong(sample, info_obj, filtered_folder, genome_size, flag):
        # I/O
        if info_obj.nanopore.trimmed:
            input_fastq = info_obj.nanopore.trimmed
        else:
            input_fastq = info_obj.nanopore.raw

        filtered_fastq = filtered_folder + sample + '.fastq.gz'

        if not os.path.exists(flag):
            cmd = ['filtlong',
                   '--keep_percent', str(95),  # drop bottom 5% reads
                   '--min_length', str(500),  # remove rejected reads from targeted sequencing
                   input_fastq]
            if genome_size:
                cmd += ['--target_bases', str(genome_size * 100)]  # keep top 100X if more reads

            # Filtlong writes to stdout
            print('\t{}'.format(sample))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            with gzip.open(filtered_fastq, 'wb') as f:
                f.write(p.communicate()[0])

        return sample, filtered_fastq

    @staticmethod
    def run_filtlong_parallel(sample_dict, output_folder, genome_size, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, genome_size, flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.run_filtlong(*x), args):
                sample_dict[results[0]].nanopore.filtered = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping filtering long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def fasta_length(input_fasta):
        with gzip.open(input_fasta, 'rt') if input_fasta.endswith('.gz') else open(input_fasta, 'r') as f:
            # Create iterator in case there are more than 1 contig in reference genome
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == '>'))

            total_len = 0
            for header in faiter:
                # Join all sequence lines of fasta entry in one line, measure its length
                # and add it to length of any other previous sequences if present
                total_len += len(''.join(s.rstrip() for s in faiter.__next__()))

            return total_len

    @staticmethod
    def assemble_flye(sample, info_obj, output_folder, genome_size, min_size, cpu, flag):
        # I/O
        assemblies_folder = output_folder + 'all_assemblies/'
        output_assembly = assemblies_folder + sample + '.fasta'

        # Create a subfolder for each sample
        output_subfolder = output_folder + sample + '/'

        # Figure out which reads we need to use
        if info_obj.nanopore.filtered:
            input_fastq = info_obj.nanopore.filtered
        elif info_obj.nanopore.trimmed:
            input_fastq = info_obj.nanopore.trimmed
        else:
            input_fastq = info_obj.nanopore.raw

        if not os.path.exists(flag):
            cmd_flye = ['flye',
                        '--nano-hq', input_fastq,
                        '--threads', str(cpu),
                        '--out-dir', output_subfolder,
                        '--iterations', str(3)]
            if min_size:
                cmd_flye += ['--min-overlap', str(min_size)]
            if genome_size:
                cmd_flye += ['--genome-size', str(genome_size)]

            print('\t{}'.format(sample))
            subprocess.run(cmd_flye, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # Rename and move assembly file
            Methods.make_folder(assemblies_folder)
            Methods.make_folder(output_subfolder)
            if os.path.exists(output_subfolder + 'assembly.fasta'):
                move(output_subfolder + 'assembly.fasta', output_assembly)

                # Rename and move assembly graph file
                assembly_graph_folder = output_folder + 'assembly_graphs/'
                Methods.make_folder(assembly_graph_folder)
                move(output_subfolder + 'assembly_graph.gfa', assembly_graph_folder + sample + '_graph.gfa')

                # Assembly graph
                cmd_bandage = ['Bandage', 'image',
                               '{}'.format(assembly_graph_folder + sample + '_graph.gfa'),
                               '{}'.format(assembly_graph_folder + sample + '_graph.png')]
                subprocess.run(cmd_bandage)
            else:
                warnings.warn('No assembly for {}!'.format(sample))

        return sample, output_assembly

    @staticmethod
    def assemble_flye_parallel(sample_dict, output_folder, genome_size, min_size, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, genome_size, min_size, int(cpu / parallel), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.assemble_flye(*x), args):
                sample_dict[results[0]].assembly = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping filtering long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def flye_assembly_stats(assembly_folder, output_folder):
        # Output file
        output_stats_file = output_folder + '/flye_stats.tsv'

        # Pandas data frame to save values
        df = pd.DataFrame(columns=['Sample', 'ReadsTotalBases', 'ReadsN50', 'AssemblyLength', 'Contigs', 'Coverage'])

        # Find log file(s) and parse info of interest
        log_list = glob(assembly_folder + '/**/flye.log', recursive=True)

        for log_file in log_list:
            sample = log_file.split('/')[-2]
            with open(log_file, 'r') as f:
                read_len = 0
                n50 = 0
                assembly_len = 0
                contigs = 0
                mean_cov = 0

                for line in f:
                    line = line.rstrip()
                    if 'Total read length:' in line:
                        read_len = line.split()[-1]
                    elif 'N50/N90' in line:
                        n50 = line.split(':')[-1].split('/')[0].strip()
                    elif 'Total length:' in line:
                        assembly_len = line.split('\t')[-1]
                    elif 'Fragments:' in line:
                        contigs = line.split('\t')[-1]
                    elif 'Mean coverage:' in line:
                        mean_cov = line.split('\t')[-1]

                data_dict = {'Sample': [sample],
                             'TotalBases': [read_len],
                             'ReadsN50': [n50],
                             'AssemblyLength': [assembly_len],
                             'Contigs': [contigs],
                             'Coverage': [mean_cov]}

                df = pd.concat([df, pd.DataFrame.from_dict(data_dict)], axis='index', ignore_index=True)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

    @staticmethod
    def assemble_shasta(sample, info_obj, output_folder, min_size, cpu, flag):

        # I/O
        # Figure out which reads we need to use
        if info_obj.nanopore.filtered:
            input_fastq = info_obj.nanopore.filtered
        elif info_obj.nanopore.trimmed:
            input_fastq = info_obj.nanopore.trimmed
        else:
            input_fastq = info_obj.nanopore.raw

        # Unzipped fastq file (needed for shasta)
        unzipped_fastq = output_folder + sample + '.fastq'

        assemblies_folder = output_folder + 'all_assemblies/'
        output_assembly = assemblies_folder + sample + '.fasta'
        output_subfolder = output_folder + sample + '/'  # Create a subfolder for each sample

        if not os.path.exists(flag):
            cmd_ungzip = ['pigz', '-dkc', input_fastq]  # To stdout
            cmd_shasta = ['shasta',
                          '--config', 'Nanopore-Oct2021',
                          '--input', unzipped_fastq,
                          '--assemblyDirectory', output_subfolder,
                          '--command', 'assemble',
                          '--threads', str(cpu)]
            if min_size:
                cmd_shasta += ['--Reads.minReadLength', str(min_size)]

            cmd_shasta_clean = ['shasta',
                                '--assemblyDirectory', output_subfolder,
                                '--command', 'cleanupBinaryData']

            # Decompress fastq for shasta
            print('\t{}'.format(sample))
            with open(unzipped_fastq, 'w') as f:
                subprocess.run(cmd_ungzip, stdout=f)

            # Run shasta assembler
            shasta_stdout = output_folder + sample + '_shasta_stdout.txt'
            with open(shasta_stdout, 'w') as f:
                subprocess.run(cmd_shasta, stdout=f, stderr=subprocess.DEVNULL)

            # Cleanup temporary files
            subprocess.run(cmd_shasta_clean, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            os.remove(unzipped_fastq)

            # Rename and move assembly file
            Methods.make_folder(assemblies_folder)
            if os.path.exists(output_subfolder + 'Assembly.fasta'):
                move(output_subfolder + 'Assembly.fasta', output_assembly)

                # Rename and move assembly graph file
                assembly_graph_folder = output_folder + 'assembly_graphs/'
                Methods.make_folder(assembly_graph_folder)
                move(output_subfolder + 'Assembly.gfa', assembly_graph_folder + sample + '.gfa')

                # Assembly graph
                cmd_bandage = ['Bandage', 'image',
                               '{}'.format(assembly_graph_folder + sample + '.gfa'),
                               '{}'.format(assembly_graph_folder + sample + '.png')]
                subprocess.run(cmd_bandage)
            else:
                warnings.warn('No assembly for {}!'.format(sample))

        return sample, output_assembly

    @staticmethod
    def assemble_shasta_parallel(sample_dict, output_folder, min_size, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, min_size, int(cpu / parallel), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.assemble_shasta(*x), args):
                sample_dict[results[0]].assembly = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping filtering long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def shasta_assembly_stats(assembly_folder, output_folder):
        # Output file
        output_stats_file = output_folder + '/shasta_stats.tsv'

        # Pandas data frame to save values
        df = pd.DataFrame(columns=['Sample', 'ReadNumber', 'TotalBases', 'ReadsN50', 'AssemblyLength', 'Contigs'])

        # Find log file(s) and parse info of interest
        log_list = glob(assembly_folder + '/*_shasta_stdout.txt', recursive=False)

        for log_file in log_list:
            sample = os.path.basename(log_file).replace("_shasta_stdout.txt", '')
            with open(log_file, 'r') as f:
                read_cnt = 0
                read_len = 0
                n50 = 0
                assembly_len = 0
                contigs = 0

                for line in f:
                    line = line.rstrip()

                    if 'Total number of reads is' in line:
                        read_cnt = line.split()[-1].replace('.', '')
                    elif 'Total number of raw bases is' in line:
                        read_len = line.split()[-1].replace('.', '')
                    elif 'N50 for read length is' in line:
                        n50 = line.split()[-2]
                    elif 'Total length of assembled sequence is' in line:
                        assembly_len = line.split()[-1]
                    elif 'The assembly graph has' in line:
                        contigs = line.split()[-3]

                data_dict = {'Sample': [sample],
                             'ReadNumber': [read_cnt],
                             'TotalBases': [read_len],
                             'ReadsN50': [n50],
                             'AssemblyLength': [assembly_len],
                             'Contigs': [contigs]}

                df = pd.concat([df, pd.DataFrame.from_dict(data_dict)], axis='index', ignore_index=True)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

        # Remove log files
        for log_file in log_list:
            os.remove(log_file)

    @staticmethod
    def run_minimap2(sample, ref, fastq_file, cpu, output_folder, keep_bam):
        print('\t{}'.format(sample))

        output_bam = output_folder + sample + '.bam'

        minimap2_cmd = ['minimap2',
                        '-a',
                        '-x', 'map-ont',
                        '-t', str(cpu),
                        '--MD',
                        '--secondary=no',
                        ref,
                        fastq_file]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '--reference', ref,
                             '-']
        samtools_markdup_cmd = ['samtools', 'markdup',
                                '-r',
                                '-@', str(cpu),
                                '-',
                                output_bam]
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index',
                              output_bam]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p4.communicate()

        # Index bam file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_minimap2_parallel(output_folder, ref, sample_dict, cpu, parallel, keep_bam):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder, keep_bam)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.run_minimap2(*x), args):
                pass