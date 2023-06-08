import os
import io
import glob
import stat
import shutil
import subprocess
from concurrent import futures
from bga_methods import Methods
import csv
import plotly.offline as offline
import plotly.graph_objs as go
import gzip
from itertools import groupby


class AssemblyQcMethods(object):
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
            for results in executor.map(lambda x: AssemblyQcMethods.run_minimap2(*x), args):
                pass

    @staticmethod
    def run_last(sample, ref_folder, query_folder, output_folder, cpu):
        # I/O
        ref_file = ref_folder + 'all_assemblies/' + sample + '.fasta'
        ref_name = '.'.join(os.path.basename(ref_file).split('.')[:-1])
        last_db = output_folder + ref_name + '.lastdb'
        query_file = query_folder + sample + '/' + sample + '.fasta'
        last_alignment = output_folder + sample + '.maf'
        dot_plot_file = output_folder + sample + '_' + os.path.dirname(ref_file).split('/')[-2] \
                        + '_vs_' + os.path.dirname(query_file).split('/')[-2] + '.png'

        # https://home.cc.umanitoba.ca/~psgendb/tutorials/bioLegato/getgenome/Last.html
        cmd_lastdb = ['lastdb', '-P', str(cpu), last_db, ref_file]
        cmd_lastal = ['lastal', '-P', str(cpu), last_db, query_file]
        cmd_dotplot = ['last-dotplot', last_alignment, dot_plot_file]

        subprocess.run(cmd_lastdb)  # Create a DB with the ref
        p = subprocess.Popen(cmd_lastal, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(last_alignment, 'w') as f:
            f.write(p.communicate()[0].decode('utf-8'))
        subprocess.run(cmd_lastal)  # Compare the polished (query) to the raw (ref)
        subprocess.run(cmd_dotplot)  # Create dotplot

    @staticmethod
    def run_last_parallel(sample_dict, ref_folder, query_folder, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref_folder, query_folder, output_folder, int(cpu / parallel))
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: AssemblyQcMethods.run_last(*x), args):
                pass

    @staticmethod
    def run_nucmer_medaka(sample, info_obj, output_folder, cpu):
        # https://mummer4.github.io/tutorial/tutorial.html
        # https://pypi.org/project/mummer-idotplot/

        # I/O
        ref = info_obj.assembly.raw
        if info_obj.assembly.medaka:
            query = info_obj.assembly.medaka
        else:  # no assembly
            return  # skip

        ouput_name = output_folder + sample + '_medaka'

        # cmd_mummer = ['mummer', '-mum',
        #               '-threads', str(cpu),
        #               '-qthreads', str(cpu),
        #               '-b', '-c', ref, query]
        # mummer_out = output_folder + sample + '.mums'
        # p = subprocess.Popen(cmd_mummer, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        # with open(mummer_out, 'w') as f:
        #     f.write(p.communicate()[0].decode('utf-8'))

        # To plot with plotly
        # cmd_mummer = ['mummer', '-maxmatch',
        #               '-F', '-L', '-b', '-l', str(10),
        #               '-threads', str(cpu),
        #               '-qthreads', str(cpu),
        #               ref, query]
        #
        # mummer_out = output_folder + sample + '.mums'
        # p = subprocess.Popen(cmd_mummer, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        # with open(mummer_out, 'w') as f:
        #     f.write(p.communicate()[0].decode('utf-8'))

        # To plot with GNUplot
        cmd_nucmer = ['nucmer',
                      '-t', str(cpu),
                      '-p', ouput_name,  # Will append '.delta'
                      '--mincluster={}'.format(100),
                      ref, query]

        subprocess.run(cmd_nucmer)

        cmd_plot = ['mummerplot',
                    '-x', '[0,{}]'.format(AssemblyQcMethods.fasta_length(ref)),
                    '-y', '[0,{}]'.format(AssemblyQcMethods.fasta_length(query)),
                    '-postscript',
                    '-p', ouput_name,
                    ouput_name + '.delta']  # nucmer output
                    # mummer_out]  # mummer output

        # GNUplot path is not defined in mummerplot when installed via conda
        Methods.fix_mummerplot()

        # Make plot
        subprocess.run(cmd_plot, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Cleanup files
        ext = ['.delta', '.fplot', '.gp', '.rplot']
        for i in ext:
            for j in glob.glob(output_folder + '/*' + i):
                if os.path.exists(j):
                    os.remove(j)

    @staticmethod
    def run_nucmer_medaka_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, int(cpu / parallel))
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: AssemblyQcMethods.run_nucmer_medaka(*x), args):
                pass

    @staticmethod
    def run_nucmer_polypolish(sample, info_obj, output_folder, cpu):
        # https://mummer4.github.io/tutorial/tutorial.html
        # https://pypi.org/project/mummer-idotplot/

        # I/O
        ref = info_obj.assembly.medaka
        query = info_obj.assembly.polypolish
        ouput_name = output_folder + sample + '_medaka'

        cmd_nucmer = ['nucmer',
                      '-t', str(cpu),
                      '-p', ouput_name,  # Will append '.delta'
                      '--mincluster={}'.format(100),
                      ref, query]
        subprocess.run(cmd_nucmer)

        cmd_plot = ['mummerplot',
                    '-x', '[0,{}]'.format(AssemblyQcMethods.fasta_length(ref)),
                    '-y', '[0,{}]'.format(AssemblyQcMethods.fasta_length(query)),
                    '-postscript',
                    '-p', ouput_name,
                    ouput_name + '.delta']  # nucmer output

        # GNUplot path is not defined in mummerplot when installed via conda
        Methods.fix_mummerplot()

        # Make plot
        subprocess.run(cmd_plot, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Cleanup files
        ext = ['.delta', '.fplot', '.gp', '.rplot']
        for i in ext:
            for j in glob.glob(output_folder + '/*' + i):
                if os.path.exists(j):
                    os.remove(j)

    @staticmethod
    def run_nucmer_polypolish_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, info_obj, output_folder, int(cpu / parallel))
                    for sample, info_obj in sample_dict.items())
            for results in executor.map(lambda x: AssemblyQcMethods.run_nucmer_polypolish(*x), args):
                pass

    @staticmethod
    def short_read_coverage(genome, r1, r2, polished_folder, cpu, sample):
        # I/O
        output_bam = polished_folder + sample + '.bam'

        # Index reference genome
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        bwa_cmd = ['bwa', 'mem', '-t', str(cpu), genome, r1, r2]
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-F', '4', '-h', '-T', genome, '-']
        samtools_fixmate_cmd = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '--reference', genome, '-']
        samtools_markdup_cmd = ['samtools', 'markdup', '-r', '-@', str(cpu), '-', output_bam]
        samtools_index_cmd = ['samtools', 'index', output_bam]  # samtools can only index chromosomes up to 512M bp.

        p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_fixmate_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p5 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p4.stdout.close()
        p5.communicate()

        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Get coverage
        cov_file = polished_folder + 'short_read_cov.tsv'
        cmd_cov = ['samtools', 'depth', output_bam]
        p = subprocess.Popen(cmd_cov, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        cov_list = list()
        for line in io.TextIOWrapper(p.stdout, encoding="utf-8"):
            cov_list.append(line.split('\t')[2])

        with open(cov_file, 'wa') as f:
            avg_cov = sum(cov_list) / len(cov_list)
            f.write('{}\t{}\n'.format(sample, round(avg_cov)))

        # Cleanup
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.bam']
        for i in ext:
            os.remove(polished_folder + "*" + i)

    @staticmethod
    def run_quast():
        pass
