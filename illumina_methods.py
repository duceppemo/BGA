import os
import io
import glob
import shutil
import subprocess
from concurrent import futures
from bga_methods import Methods


class IlluminaMethods(object):
    @staticmethod
    def trim_illumina_fastp_paired(sample, info_obj, illumina_trimmed_folder, report_folder, cpu, flag):
        # I/O
        try:
            in_r1 = info_obj.illumina.raw[0]
            in_r2 = info_obj.illumina.raw[1]

            out_r1 = illumina_trimmed_folder + sample + '_R1.fastq.gz'
            out_r2 = illumina_trimmed_folder + sample + '_R2.fastq.gz'

            if not os.path.exists(flag):
                # Check that we have paired-end reads
                if os.path.exists(in_r1) and os.path.exists(in_r2):
                    cmd = ['fastp',
                           '--in1', in_r1,
                           '--in2', in_r2,
                           '--out1', out_r1,
                           '--out2', out_r2,
                           '--detect_adapter_for_pe',
                           '--cut_right',
                           '--cut_right_mean_quality', str(10),
                           '--length_required', str(64),
                           '--html', report_folder + sample + '.html',
                           '--thread', str(cpu)]

                    print('\t{}'.format(sample))
                    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                else:
                    print('Paired-end data missing for {}. Skipping short read trimming.'.format(sample))
                    out_r1 = ''
                    out_r2 = ''
        except IndexError:
            out_r1 = ''
            out_r2 = ''

        return sample, out_r1, out_r2

    @staticmethod
    def trim_illumina_fastp_paired_parallel(sample_dict, output_folder, report_folder, cpu, parallel, flag):
        Methods.make_folder(output_folder)
        Methods.make_folder(report_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, report_folder,
                     int(cpu / parallel), flag) for sample, sample_obj in sample_dict.items())
            for results in executor.map(lambda x: IlluminaMethods.trim_illumina_fastp_paired(*x), args):
                sample_dict[results[0]].illumina.trimmed.insert(0, results[1])
                sample_dict[results[0]].illumina.trimmed.insert(0, results[2])

        if os.path.exists(flag):  # Already performed
            print('\tSkipping trimming short reads. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)
            os.remove('fastp.json')

    @staticmethod
    def map_minimap2_short_pe(genome, r1, r2, output_folder, cpu, sample):
        # I/O
        out_bam = output_folder + sample + '.bam'

        minimap2_cmd = ['minimap2', '-a', '-x', 'sr', '-t', str(cpu), genome, r1, r2]
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-F', '4', '-h', '-']
        samtools_fixmate_cmd = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '-']
        samtools_markdup_cmd = ['samtools', 'markdup', '-r', '-@', str(cpu), '-', out_bam]
        samtools_index_cmd = ['samtools', 'index', out_bam]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_fixmate_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_sort_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p5 = subprocess.Popen(samtools_markdup_cmd, stdin=p4.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p4.stdout.close()
        p5.communicate()

        # Index bam file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_bwa_paired(genome, r1, r2, output_folder, cpu, sample):
        # I/O
        out_bam = output_folder + sample + '.bam'

        # Index reference genome
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        # Map reads
        cmd_bwa_mem = ['bwa', 'mem', '-t', str(cpu), genome, r1, r2]
        cmd_samtools_view = ['samtools', 'view', '-@', str(cpu), '-F', str('0x4'), '-b', '-']  # only keep mapped
        cmd_samtools_fixmate = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        cmd_samtools_sort = ['samtools', 'sort', '-@', str(cpu), '-']
        cmd_samtools_markdup = ['samtools', 'markdup', '-@', str(cpu), '-r', '-', out_bam]

        p1 = subprocess.Popen(cmd_bwa_mem, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(cmd_samtools_view, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(cmd_samtools_fixmate, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(cmd_samtools_sort, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p5 = subprocess.Popen(cmd_samtools_markdup, stdin=p4.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p4.stdout.close()
        p5.communicate()

        # Index bam file
        cmd = ['samtools', 'index', '-@', str(cpu), out_bam]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def fix_fasta(input_fasta, output_fasta):
        cmd = ['seqtk', 'seq',
               input_fasta]

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        with open(output_fasta, 'w') as f:
            for line in io.TextIOWrapper(p.stdout, encoding="utf-8"):
                if line.startswith('>'):
                    f.write(line)
                else:
                    f.write(line.upper())

    @staticmethod
    # def run_polypolish(genome, r1, r2, polished_folder, sample, n_run, cpu):
    def run_polypolish(genome, r1, r2, polished_folder, sample, cpu):
        # I/O
        sam_r1 = polished_folder + sample + '_R1.sam'
        sam_r2 = polished_folder + sample + '_R2.sam'

        # Index genome
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Map reads
        cmd_map_r1 = ['bwa', 'mem', '-t', str(cpu), '-a', genome, r1]
        cmd_map_r2 = ['bwa', 'mem', '-t', str(cpu), '-a', genome, r2]

        p1 = subprocess.Popen(cmd_map_r1, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r1, 'w') as f:
            f.write(p1.communicate()[0].decode('utf-8'))

        p2 = subprocess.Popen(cmd_map_r2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r2, 'w') as f:
            f.write(p2.communicate()[0].decode('utf-8'))

        cmd_pp = ['polypolish', 'polish', genome, sam_r1, sam_r2]
        p3 = subprocess.Popen(cmd_pp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sdtout, stderr = p3.communicate()

        fasta_pp = polished_folder + sample + '.polypolish.fasta'
        # polypolish_log = polished_folder + sample + '.polypolish' + '_round' + str(n_run + 1) + '.log'
        polypolish_log = polished_folder + sample + '.polypolish.log'
        with open(fasta_pp, 'w') as f:
            f.write(sdtout.decode('utf-8'))
        with open(polypolish_log, 'w') as f:
            f.write(stderr.decode('utf-8'))

        # Fix fasta
        fixed_fasta = polished_folder + sample + '.fasta'
        IlluminaMethods.fix_fasta(fasta_pp, fixed_fasta)

        # Cleanup
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.sam']
        for i in ext:
            for j in glob.glob(polished_folder + '/' + sample + '*' + i):
                if os.path.exists(j):
                    os.remove(j)
        os.remove(fasta_pp)

        return fixed_fasta

    @staticmethod
    def run_pypolca(genome, r1, r2, polished_folder, sample, cpu):
        # https://github.com/gbouras13/pypolca
        # The polished output FASTA will be {prefix}_corrected.fasta in the specified output directory
        # and the POLCA report will be the textfile {prefix}.report
        tmp_folder = polished_folder + sample + '/'
        fasta_pp = tmp_folder + sample + '_corrected.fasta'
        Methods.make_folder(tmp_folder)

        cmd = ['pypolca', 'run',
               '--assembly', genome,
               '--reads1', r1,
               '--reads2', r2,
               '--output', tmp_folder,
               '--prefix', sample,
               '--careful',
               '--force',
               '--threads', str(cpu)]

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()

        pypolca_log = polished_folder + sample + '.pypolca.log'
        with open(pypolca_log, 'w') as f:
            f.write(stdout.decode('utf-8'))

        # Fix fasta
        fixed_fasta = polished_folder + sample + '.fasta'
        IlluminaMethods.fix_fasta(fasta_pp, fixed_fasta)

        # Cleanup
        try:
            # If already ran
            if os.path.exists(polished_folder + sample + '.report'):
                os.remove(polished_folder + sample + '.report')
            shutil.move(tmp_folder + sample + '.report', polished_folder)
            if os.path.exists(polished_folder + sample + '.vcf'):
                os.remove(polished_folder + sample + '.vcf')
            shutil.move(tmp_folder + sample + '.vcf', polished_folder)  # Sometimes not present
        except FileNotFoundError:
            pass
        shutil.rmtree(tmp_folder)

        return fixed_fasta

    @staticmethod
    def fix_fasta_header(genome):
        # Rename fasta file
        tmp_file = genome + '.tmp'
        os.rename(genome, genome + '.tmp')

        with open(genome, 'w') as out_fh:
            with open(tmp_file, 'r') as in_fh:
                i = 0
                for line in in_fh:
                    if line.startswith('>'):
                        i += 1
                        out_fh.write('>contig_{}\n'.format(i))
                    else:
                        out_fh.write(line)

        # Cleanup temp file
        os.remove(tmp_file)

    @staticmethod
    def polish(sample, info_obj, output_folder, cpu, flag):
        # I/O
        try:
            r1 = info_obj.illumina.trimmed[0]
            r2 = info_obj.illumina.trimmed[1]
        except IndexError:
            r1 = info_obj.illumina.raw[0]
            r2 = info_obj.illumina.raw[1]

        genome = info_obj.assembly.medaka
        polished_assembly = output_folder + sample + '.fasta'

        if not os.path.exists(flag):
            if os.path.exists(genome):
                if os.path.exists(r1) and os.path.exists(r2):
                    # Check if there is an assembly available for that sample
                    print('\t{}'.format(sample))

                    # Create polish folder
                    Methods.make_folder(output_folder)

                    # Polypolish
                    genome = IlluminaMethods.run_polypolish(genome, r1, r2, output_folder, sample, cpu)

                    # Pypolca
                    IlluminaMethods.run_pypolca(genome, r1, r2, output_folder, sample, cpu)
                else:
                    print('Paired-end data missing for {}. Skipping short read polishing.'.format(sample))
                    polished_assembly = ''
            else:
                print('No assembly for {}. Skipping short read polishing.'.format(sample))
                polished_assembly = ''

        # Need this in case a file is missing and the pipeline is skipping already completed steps
        if os.path.exists(polished_assembly):
            IlluminaMethods.fix_fasta_header(polished_assembly)

        return sample, polished_assembly

    @staticmethod
    def polish_parallel(sample_dict, output_folder, cpu, parallel, flag):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, output_folder, int(cpu / parallel), flag)
                    for sample, sample_obj in sample_dict.items())
            for results in executor.map(lambda x: IlluminaMethods.polish(*x), args):
                sample_dict[results[0]].assembly.polypolish = results[1]

        if os.path.exists(flag):  # Already performed
            print('\tSkipping short read polishing. Already done.')
        else:  # Create the done flag
            Methods.flag_done(flag)
