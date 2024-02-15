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
    def run_nextpolish(genome, r1, r2, polished_folder, cpu, sample):
        # Index genome
        subprocess.run(['samtools', 'faidx', genome], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Map paired-end reads to genome
        IlluminaMethods.map_bwa_paired(genome, r1, r2, polished_folder, cpu, sample)
        # print('\t\tMapping reads')
        # IlluminaMethods.map_minimap2_short_pe(genome, r1, r2, polished_folder, cpu, sample)

        # nextPolish subscripts location
        # ./bin/nextPolish
        # ./share/nextpolish-1.4.1/lib/nextpolish1.py
        np_path = shutil.which('nextPolish')
        np1_path = '/'.join(np_path.split('/')[:-2]) + '/share/nextpolish-1.4.1/lib/nextpolish1.py'

        # Polish
        cmd = ['python', np1_path,
               '--genome', genome,
               '--task', str(1),
               '--process', str(cpu),
               '--uppercase',
               '--bam_sgs', polished_folder + sample + '.bam',
               '-ploidy', str(1)]

        fasta_np = polished_folder + sample + '.nextpolish.fasta'
        fasta_fixed = polished_folder + sample + '.fasta'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)  # Debug
        # print('\t\tRunning NextPolish')
        with open(fasta_np, 'w') as f:
            f.write(p.communicate()[0].decode('utf-8'))

        # Fix fasta
        # print('\t\tFixing fasta')
        IlluminaMethods.fix_fasta(fasta_np, fasta_fixed)

        # Cleanup
        genome_folder = os.path.dirname(genome)
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.fai', '.bam', '.bai']
        for i in ext:
            for j in glob.glob(genome_folder + '/*' + i):
                if os.path.exists(j):
                    os.remove(j)
        os.remove(fasta_np)

        return fasta_fixed

    @staticmethod
    def run_ntedit(genome, r1, r2, polished_folder, cpu, sample):
        cmd_nthits = ['nthits',
                      '-t', str(cpu),
                      '-k', str(40),
                      '-p', polished_folder + sample,
                      '--outbloom',
                      '--solid',
                      r1, r2]

        cmd_ntedit = ['ntedit',
                      '-t', str(cpu),
                      '-f', genome,
                      '-b', polished_folder + sample,
                      '-r', polished_folder + sample + '_k40.bf',
                      '-m', str(1)]

        subprocess.run(cmd_nthits, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.run(cmd_ntedit, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Fix fasta
        fixed_fasta = polished_folder + sample + '.fasta'
        IlluminaMethods.fix_fasta(polished_folder + sample + '_edited.fa', fixed_fasta)

        # Cleanup
        ext = ['_variants.vcf', '_changes.tsv', '_k40.bf', '_edited.fa']
        for i in ext:
            for j in glob.glob(polished_folder + '/*' + i):
                if os.path.exists(j):
                    os.remove(j)

        return fixed_fasta

    @staticmethod
    def run_polypolish(genome, r1, r2, polished_folder, sample):
        # I/O
        sam_r1 = polished_folder + sample + '_R1.sam'
        sam_r2 = polished_folder + sample + '_R2.sam'

        # Index genome
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Map reads
        cmd_map_r1 = ['bwa', 'mem', '-a', genome, r1]
        cmd_map_r2 = ['bwa', 'mem', '-a', genome, r2]

        p1 = subprocess.Popen(cmd_map_r1, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r1, 'w') as f:
            f.write(p1.communicate()[0].decode('utf-8'))

        p2 = subprocess.Popen(cmd_map_r2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r2, 'w') as f:
            f.write(p2.communicate()[0].decode('utf-8'))

        cmd_pp = ['polypolish', 'polish', genome, sam_r1, sam_r2]
        fasta_pp = polished_folder + sample + '.polypolish.fasta'
        p3 = subprocess.Popen(cmd_pp, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(fasta_pp, 'w') as f:
            f.write(p3.communicate()[0].decode('utf-8'))

        # Fix fasta
        fixed_fasta = polished_folder + sample + '.fasta'
        IlluminaMethods.fix_fasta(fasta_pp, fixed_fasta)

        # Cleanup
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.sam']
        for i in ext:
            for j in glob.glob(polished_folder + '/*' + i):
                if os.path.exists(j):
                    os.remove(j)
        os.remove(fasta_pp)

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

        # TODO: add flag check
        if not os.path.exists(flag):
            if os.path.exists(genome):
                if os.path.exists(r1) and os.path.exists(r2):
                    # Check if there is an assembly available for that sample
                    print('\t{}'.format(sample))
                    # Create polish folder
                    Methods.make_folder(output_folder)

                    # NextPolish
                    # print('\tNextPolish 1st round')  # Debug
                    genome = IlluminaMethods.run_nextpolish(genome, r1, r2, output_folder, cpu, sample)
                    # print('\tNextPolish 2nd round')  # Debug
                    genome = IlluminaMethods.run_nextpolish(genome, r1, r2, output_folder, cpu, sample)

                    # ntEdit
                    # print('\tntEdit')  # Debug
                    genome = IlluminaMethods.run_ntedit(genome, r1, r2, output_folder, cpu, sample)

                    # Polypolish
                    for i in range(3):
                        # print('\tPolypolish round {}'.format(i))  # Debug
                        genome = IlluminaMethods.run_polypolish(genome, r1, r2, output_folder, sample)
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
