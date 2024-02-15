import os
import gzip
import shutil
import subprocess
import pandas as pd
from glob import glob
from shutil import move
from concurrent import futures
from bga_methods import Methods, Sample, Assembly


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
            print('\tSkipping trimming long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def run_filtlong(sample, info_obj, filtered_folder, genome_size, flag):
        # I/O
        try:
            input_fastq = info_obj.nanopore.trimmed
        except AttributeError:
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

        # Need this in case a file is missing and the pipeline is skipping already completed steps
        if not os.path.exists(filtered_fastq):
            filtered_fastq = ''

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
    def assemble_flye(sample, info_obj, output_folder, gfa_folder, genome_size, min_size, cpu, flag):
        # I/O
        output_assembly = output_folder + sample + '.fasta'
        output_subfolder = output_folder + sample + '/'  # Create a subfolder for each sample

        # Figure out which reads we need to use
        try:
            input_fastq = info_obj.nanopore.filtered
        except AttributeError:
            try:
                input_fastq = info_obj.nanopore.trimmed
            except AttributeError:
                input_fastq = info_obj.nanopore.raw

        if not os.path.exists(flag):  # Skip if already preformed (flag file present)
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
            Methods.make_folder(output_subfolder)
            subprocess.run(cmd_flye, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # Rename and move assembly file
            if os.path.exists(output_subfolder + 'assembly.fasta'):  # Assembly was successful
                move(output_subfolder + 'assembly.fasta', output_assembly)

                # Rename and move assembly graph file
                Methods.make_folder(gfa_folder)
                move(output_subfolder + 'assembly_graph.gfa', gfa_folder + sample + '.gfa')

                # Assembly graph
                cmd_bandage = ['Bandage', 'image',
                               '{}'.format(gfa_folder + sample + '.gfa'),
                               '{}'.format(gfa_folder + sample + '.png')]
                subprocess.run(cmd_bandage, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            else:  # No assembly output
                output_assembly = ''

        # Need this in case a file is missing and the pipeline is skipping already completed steps
        if not os.path.exists(output_assembly):
            output_assembly = ''

        return sample, output_assembly

    @staticmethod
    def assemble_flye_parallel(sample_dict, output_folder, gfa_folder, genome_size, min_size, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, gfa_folder, genome_size, min_size, int(cpu / parallel), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.assemble_flye(*x), args):
                # Initiate object
                my_sample = Sample()
                my_sample.assembly = Assembly()

                sample_dict[results[0]].assembly.raw = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping assembling long reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)

    @staticmethod
    def flye_assembly_stats(assembly_folder, output_folder):
        # Output file
        output_stats_file = output_folder + '/flye_stats.tsv'

        # Pandas data frame to save values
        df = pd.DataFrame(columns=['Sample', 'TotalReadLength', 'ReadsN50', 'AssemblyLength', 'Contigs', 'Coverage'])

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
                             'TotalReadLength': [read_len],
                             'ReadsN50': [n50],
                             'AssemblyLength': [assembly_len],
                             'Contigs': [contigs],
                             'Coverage': [mean_cov]}

                df = pd.concat([df, pd.DataFrame.from_dict(data_dict)], axis='index', ignore_index=True)

            # Remove assembly folder
            sub_folder = os.path.dirname(log_file)
            shutil.rmtree(sub_folder)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

    @staticmethod
    def assemble_shasta(sample, info_obj, output_folder, gfa_folder, min_size, cpu, flag):
        # I/O
        # Figure out which reads we need to use
        try:
            input_fastq = info_obj.nanopore.filtered
        except AttributeError:
            try:
                input_fastq = info_obj.nanopore.trimmed
            except AttributeError:
                input_fastq = info_obj.nanopore.raw

        # Unzipped fastq file (needed for shasta)
        unzipped_fastq = output_folder + sample + '.fastq'
        fasta = output_folder + sample + '.fasta'

        output_assembly = output_folder + sample + '.fasta'
        output_subfolder = output_folder + sample + '/'  # Create a subfolder for each sample

        if not os.path.exists(flag):
            cmd_ungzip = ['pigz', '-dkc', input_fastq]  # To stdout

            cmd_fastq_to_fasta = ['seqtk', 'seq',
                                  '-a', input_fastq]

            cmd_shasta = ['shasta',
                          '--config', 'Nanopore-May2022',
                          # '--input', unzipped_fastq,
                          '--input', fasta,
                          '--assemblyDirectory', output_subfolder,
                          '--threads', str(cpu)]
            if min_size:
                cmd_shasta += ['--Reads.minReadLength', str(min_size)]
            else:
                cmd_shasta += ['--Reads.minReadLength', str(3000)]

            cmd_shasta_clean = ['shasta',
                                '--assemblyDirectory', output_subfolder,
                                '--command', 'cleanupBinaryData']

            # Decompress fastq for shasta
            # print('\t{}'.format(sample))
            # with open(unzipped_fastq, 'w') as f:
            #     subprocess.run(cmd_ungzip, stdout=f)

            # Convert fastq to fasta
            with open(fasta, 'w') as f:
                subprocess.run(cmd_fastq_to_fasta, stdout=f)

            # Run shasta assembler
            # Need this file to get the assembly stats
            shasta_stdout = output_folder + sample + '_shasta_stdout.txt'
            with open(shasta_stdout, 'w') as f:
                subprocess.run(cmd_shasta, stdout=f, stderr=subprocess.DEVNULL)
            # subprocess.run(cmd_shasta, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            # Cleanup temporary files
            subprocess.run(cmd_shasta_clean, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            # os.remove(unzipped_fastq)
            os.remove(fasta)

            # Rename and move assembly file
            # Assembly and gfa files are always created, but empty if no assemblies
            if os.stat(output_subfolder + 'Assembly.fasta').st_size != 0:  # If assembly file not empty
                move(output_subfolder + 'Assembly.fasta', output_assembly)

                # Rename and move assembly graph file
                Methods.make_folder(gfa_folder)
                move(output_subfolder + 'Assembly.gfa', gfa_folder + sample + '.gfa')

                # Assembly graph
                cmd_bandage = ['Bandage', 'image',
                               '{}'.format(gfa_folder + sample + '.gfa'),
                               '{}'.format(gfa_folder + sample + '.png')]
                subprocess.run(cmd_bandage, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            else:
                output_assembly = ''

        # Need this in case a file is missing and the pipeline is skipping already completed steps
        if not os.path.exists(output_assembly):
            output_assembly = ''

        return sample, output_assembly

    @staticmethod
    def assemble_shasta_parallel(sample_dict, output_folder, gfa_folder, min_size, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            # sample, info_obj, output_folder, gfa_folder, min_size, cpu, flag
            args = ((sample, info_obj, output_folder, gfa_folder, min_size, int(cpu / parallel), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.assemble_shasta(*x), args):
                sample_dict[results[0]].assembly.raw = results[1]

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

                # Remove log files
                os.remove(log_file)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

        # Delete sample subfolder
        for root, directories, filenames in os.walk(assembly_folder):
            for folder in directories:
                shutil.rmtree(os.path.join(root, folder))

    @staticmethod
    def polish_medaka(sample, info_obj, output_folder, model, cpu, flag):
        # I/O
        # Figure out which reads we need to use
        try:
            input_fastq = info_obj.nanopore.filtered
        except AttributeError:
            try:
                input_fastq = info_obj.nanopore.trimmed
            except AttributeError:
                input_fastq = info_obj.nanopore.raw

        input_assembly = info_obj.assembly.raw

        sample_subfolder = output_folder + sample + '/'
        default_medaka_output = sample_subfolder + 'consensus.fasta'
        polished_assembly = output_folder + sample + '.fasta'
        # medaka_stdout = sample_subfolder + sample + '_medaka_stdout.txt'

        cmd_medaka = ['medaka_consensus',
                      '-i', input_fastq,
                      '-d', input_assembly,
                      '-o', sample_subfolder,
                      '-t', str(cpu),
                      '-m', model]

        if not os.path.exists(flag):
            if os.path.exists(input_assembly):
                # Create output folders
                Methods.make_folder(output_folder)
                Methods.make_folder(sample_subfolder)

                print('\t{}'.format(sample))

                # Run medaka and save STDOUT to file
                # Had to do this to figure out why medaka was not completing sometimes
                # Turns out it's because BCFtools sometimes does not install properly with conda
                # with open(medaka_stdout, 'w') as f:
                #     subprocess.run(cmd_medaka, stdout=f, stderr=subprocess.STDOUT)
                subprocess.run(cmd_medaka, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                # Rename medaka output assembly
                if os.path.exists(default_medaka_output):
                    os.rename(default_medaka_output, polished_assembly)

                    # Reformat fasta to have 80 characters per line
                    Methods.format_fasta(polished_assembly, polished_assembly + '.tmp')
                    os.rename(polished_assembly + '.tmp', polished_assembly)

                # Cleanup
                shutil.rmtree(sample_subfolder)
                # ext = ['.gfa', '.log', '_medaka.fasta', '.bed', '.hdf', '.bam', '.bai']
                # for i in ext:
                #     for j in glob(sample_subfolder + '/*' + i):
                #         if os.path.exists(j):
                #             os.remove(j)
                try:
                    os.remove(input_assembly + '.map-ont.mmi')
                    os.remove(input_assembly + '.fai')
                except FileNotFoundError:
                    pass
            else:
                print('\tNo assembly for {}'.format(sample))
                polished_assembly = ''

        # Need this in case a file is missing and the pipeline is skipping already completed steps
        if not os.path.exists(polished_assembly):
            polished_assembly = ''

        return sample, polished_assembly

    @staticmethod
    def polish_medaka_parallel(sample_dict, output_folder, model, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        # Medaka crashes when multiple instances run in parallel
        # with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
        #     args = ((sample, info_obj, output_folder, model, int(cpu / parallel), flag)
        #             for sample, info_obj in sample_dict.items())
        #     for results in executor.map(lambda x: NanoporeMethods.polish_medaka(*x), args):
        #         sample_dict[results[0]].assembly.medaka = results[1]

        # Running one at the time.
        with futures.ThreadPoolExecutor(max_workers=1) as executor:
            args = ((sample, info_obj, output_folder, model, int(cpu), flag)
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: NanoporeMethods.polish_medaka(*x), args):
                sample_dict[results[0]].assembly.medaka = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping long read polishing. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)
