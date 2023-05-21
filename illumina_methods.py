import os
import shutil
import subprocess
from concurrent import futures
import io


class IlluminaMethods(object):
    @staticmethod
    def trim_illumina_fastp_paired(sample, info_obj, illumina_trimmed_folder, report_folder, cpu):
        in_r1 = info_obj.illumina.raq.r1
        in_r2 = info_obj.illumina.raw.r1

        out_r1 = illumina_trimmed_folder + sample + '_R1.fastq.gz'
        out_r2 = illumina_trimmed_folder + sample + '_R2.fastq.gz'

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

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        return sample, out_r1, out_r2

    @staticmethod
    def trim_illumina_fastp_paired_parallel(sample_dict, illumina_trimmed_folder, report_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, sample_obj, illumina_trimmed_folder, report_folder,
                     int(cpu / parallel)) for sample, sample_obj in sample_dict.items())
            for results in executor.map(lambda x: IlluminaMethods.trim_illumina_fastp_paired(*x), args):
                sample_dict[results[0]].illumina.trimmed.r1 = results[1]
                sample_dict[results[0]].illumina.trimmed.r2 = results[2]

    @staticmethod
    def map_bwa_paired(genome, info_obj, output_folder, cpu, sample):
        # Index reference genome
        cmd_bwa_index = ['bwa', 'index', genome]

        # I/O
        if info_obj.illumina.trimmed.r1:
            r1 = info_obj.illumina.trimmed.r1
            r2 = info_obj.illumina.trimmed.r2
        else:
            r1 = info_obj.illumina.raw.r1
            r2 = info_obj.illumina.raw.r2

        out_bam = output_folder + sample + '.bam'

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
        p4 = subprocess.Popen(cmd_samtools_sort, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p5 = subprocess.Popen(cmd_samtools_markdup, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p4.stdout.close()
        p5.communicate()

        # Index bam file
        cmd = ['samtools', 'index', '-@', str(cpu), out_bam]

    @staticmethod
    def fix_fasta(input_fasta, output_fasta):
        cmd = ['seqkt', 'seq',
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
        # Map paired-end reads to genome
        IlluminaMethods.map_bwa_paired(genome, r1, r2, polished_folder, cpu, sample)

        # nextPolish subscripts location
        # ./bin/nextPolish
        # ./share/nextpolish-1.4.1/lib/nextpolish1.py
        np_path = shutil.which('nextPolish')
        np1_path = '/'.join(np_path.split('/')[:-2]) + '/share/nextpolish-1.4.1/lib/nextpolish1.py'

        # Polish
        cmd = ['python', np1_path,
               '--genome', genome,
               '--task', str(1),
               '--precess', str(cpu),
               '--bam_sgs', polished_folder + sample + '.bam',
               '-ploidy', str(1)]

        fasta_np = polished_folder + sample + '.nextpolish.fasta'
        fasta_fixed = polished_folder + sample + '.fasta'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(fasta_np, 'w') as f:
            f.write(p.communicate()[0])

        # Fix fasta
        IlluminaMethods.fix_fasta(fasta_np, fasta_fixed)

        # Cleanup
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.fai', '.bam', '.bai']
        for i in ext:
            os.remove(polished_folder + "*" + i)
        os.remove(fasta_np)

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

        subprocess.run(cmd_nthits)
        subprocess.run(cmd_ntedit)

        # Fix fasta
        IlluminaMethods.fix_fasta(polished_folder + sample + '_edited.fa',
                                  polished_folder + sample + '.fasta')

        # Cleanup
        ext = ['_variants.vcf', '_changes.tsv', '_k40.bf', '_edited.fa']
        for i in ext:
            os.remove(polished_folder + "*" + i)

    @staticmethod
    def run_polypolish(genome, r1, r2, polished_folder, sample):
        # I/O
        sam_r1 = polished_folder + sample + '_R1.sam'
        sam_r2 = polished_folder + sample + '_R2.sam'

        # Index genome
        cmd_bwa_index = ['bwa', 'index', genome]
        subprocess.run(cmd_bwa_index)

        # Map reads
        cmd_map_r1 = ['bwa', 'mem', '-a', genome, r1]
        cmd_map_r2 = ['bwa', 'mem', '-a', genome, r2]

        p1 = subprocess.Popen(cmd_map_r1, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r1, 'w') as f:
            f.write(p1.communicate()[0])

        p2 = subprocess.Popen(cmd_map_r2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with open(sam_r2, 'w') as f:
            f.write(p2.communicate()[0])

        cmd_pp = ['polypolish', genome, sam_r1, sam_r2]
        fasta_pp = polished_folder + sample + '.polypolish.fasta'
        with open(fasta_pp, 'w') as f:
            f.write(p2.communicate()[0])

        # Fix fasta
        fasta_out = polished_folder + sample + '.fasta'
        IlluminaMethods.fix_fasta(fasta_pp, fasta_out)

        # Cleanup
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa', '.sam']
        for i in ext:
            os.remove(polished_folder + "*" + i)
        os.remove(fasta_pp)

    @staticmethod
    def fix_fasta_header(genome):
        # Rename fasta file
        tmp_file =  genome + '.tmp'
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
    def short_read_coverage(genome, r1, r2, polished_folder, cpu, sample):
        # I/O
        output_bam = polished_folder + sample + '.bam'

        # Index reference genome
        cmd_bwa_index = ['bwa', 'index', genome]

        bwa_cmd = ['bwa', 'mem', '-t', str(cpu), genome, r1, r2]
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-F', '4', '-h', '-T', genome, '-']
        samtools_fixmate_cmd = ['samtools', 'fixmate', '-@', str(cpu), '-m', '-', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '--reference', genome, '-']
        samtools_markdup_cmd = ['samtools', 'markdup', '-r', '-@', str(cpu), '-', output_bam]

        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index', output_bam]

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
    def polish(sample_dict, polished_folder, cpu):
        for sample, info_obj in sample_dict.items():
            # I/O
            if info_obj.illumina.trimmed.r1:
                r1 = info_obj.illumina.trimmed.r1
                r2 = info_obj.illumina.trimmed.r2
            else:
                r1 = info_obj.illumina.raw.r1
                r2 = info_obj.illumina.raw.r2

            # NextPolish
            for i in range(2):
                IlluminaMethods.run_nextpolish(info_obj.assembly, r1, r2, polished_folder, cpu, sample)

            # ntEdit
            IlluminaMethods.run_ntedit(info_obj.assembly, r1, r2, polished_folder, cpu, sample)

            # Polypolish
            for i in range(3):
                IlluminaMethods.run_polypolish(info_obj.assembly, r1, r2, polished_folder, sample)
