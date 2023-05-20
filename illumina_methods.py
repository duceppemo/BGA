import shutil
import subprocess
from concurrent import futures


class IlluminaMethods(object):
    @staticmethod
    def trim_illumina_fastp_paired(sample, in_r1, in_r2, illumina_trimmed_folder, report_folder, cpu):
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

    @staticmethod
    def trim_illumina_fastp_paired_parallel(sample_dict, illumina_trimmed_folder, report_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, r1, r2, illumina_trimmed_folder, report_folder, int(cpu / parallel))
                    for r1, r2 in sample_obj.illumina for sample, sample_obj in sample_dict.items())
            for results in executor.map(lambda x: IlluminaMethods.trim_illumina_fastp_paired(*x), args):
                pass

    @staticmethod
    def map_bwa_paired(genome, r1, r2, output_folder, cpu, sample):
        # Index reference genome
        cmd_bwa_index = ['bwa', 'index', genome]

        # Map reads
        out_bam = output_folder + sample + '.bam'

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
    def run_nextpolish(genome, r1, r2, output_folder, cpu, sample):
        # Map paired-end reads to genome
        IlluminaMethods.map_bwa_paired(genome, r1, r2, output_folder, cpu, sample)

        # nextPolish subscripts location
        np_path = shutil.which('nextPolish')

        # Polish
        cmd = ['nextPolish',
               '--genome', genome,
               '--task', str(1),
               '--precess', str(cpu),
               '--bam_sgs', output_folder + sample + '.bam',
               '-ploidy', str(1)]
        """
        python $HOME/prog/NextPolish/lib/nextpolish1.py \
            --genome "$genome" \
            --task 1 \
            --process $((cpu/maxProc)) \
            --bam_sgs "${out}"/"${sample}".bam \
            -ploidy 1 \
            > "${out}"/"${sample}".nextpolish.fasta

        # Make uppercase and one line per sequence
        seqtk seq "${out}"/"${sample}".nextpolish.fasta | \
            awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
            > "${out}"/"${sample}".fasta

        # Cleanup
        rm "${genome}".amb \
            "${genome}".ann \
            "${genome}".bwt \
            "${genome}".pac \
            "${genome}".sa \
            "${genome}".fai \
            "${out}"/"${sample}".bam* \
            "${out}"/"${sample}".nextpolish.fasta

        # Update genome for next step
        genome="${out}"/"${sample}".fasta
    }
        """
