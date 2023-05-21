import subprocess
import os
import sys
from concurrent import futures
import pathlib
from psutil import virtual_memory
from multiprocessing import cpu_count
import gzip
from glob import glob
import warnings


class Sample(object):
    def __getattr__(self, name):
        self.__dict__[name] = Sample()
        return self.__dict__[name]


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        # Check if exists and is a folder
        if not (os.path.exists(my_input) and os.path.isdir(my_input)):
            raise Exception('Please select an existing folder as input.')

        # List content of folder
        file_list = os.listdir(my_input)

        # Check if folder is not empty and all files have the accepted extensions
        if not all([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure all files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def check_version(log_file):
        # Not being used right now because versions are captured in the requirements.txt file
        with open(log_file, 'w') as f:
            # Python
            p = subprocess.Popen(['python', '--version'])
            stderr, stdout = p.communicate()
            # Python 3.10.11  # Outputs just this line

            # Porechop
            p = subprocess.Popen(['porechop', '--version'])
            stderr, stdout = p.communicate()
            # 0.2.4  # Outputs just this line

            # BBmap suite
            p = subprocess.Popen(['bbduk.sh', '--version'])
            stderr, stdout = p.communicate()
            # BBMap version 39.01  # Got to fetch that line from output

            # Filtlong
            p = subprocess.Popen(['filtlong', '--version'])
            stderr, stdout = p.communicate()
            # Filtlong v0.2.1  # Outputs just this line

            # SAMtools
            p = subprocess.Popen(['samtools'])
            stderr, stdout = p.communicate()
            # Version: 1.17 (using htslib 1.17)  # Got to fetch that line from output


            # Minimap2
            p = subprocess.Popen(['minimap2', '--version'])
            stderr, stdout = p.communicate()
            # 2.26-r1175  # Outputs just this line

            # bwa
            p = subprocess.Popen(['bwa'])
            stderr, stdout = p.communicate()
            # Version: 0.7.17-r1188  # Got to fetch that line from output

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_nanopore_files(in_folder):
        sample_dict = dict()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].split('_')[0]
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]

                    sample_obj = Sample()
                    sample_dict[sample] = sample_obj
                    sample_dict[sample].nanopore.raw = file_path

        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def get_illumina_files(in_folder, sample_dict):
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].split('_')[0]
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]

                    if '_R1' in filename:
                        sample_dict[sample].illumina.raw.r1 = file_path
                    elif '_R2' in filename:
                        sample_dict[sample].illumina.raw.r2 = file_path

        return sample_dict

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob(folder + '/*' + extension)

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass

    @staticmethod
    def gzipped_file_size(gzipped_file):
        with gzip.open(gzipped_file, 'rb') as f:
            return f.seek(0, whence=2)

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
        if os.path.exists(output_bam):
            if os.stat(output_bam).st_size != 0:  # bam file exists and not empty
                subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                # Convert bam to fastq
                Methods.get_fastq_from_bam(sample, output_bam, fastq_file, output_folder)

                # Remove bam
                if not keep_bam:
                    bam_list = glob(output_folder + '*.bam*')
                    for bam in bam_list:
                        os.remove(bam)
        else:
            warnings.warn('No reads were extracted for {}!'.format(sample))

    @staticmethod
    def run_minimap2_parallel(output_folder, ref, sample_dict, cpu, parallel, keep_bam):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder, keep_bam)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_minimap2(*x), args):
                pass
