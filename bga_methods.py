
import os
import sys
import gzip
import shutil
import pathlib
import textwrap
import subprocess
from glob import glob
from psutil import virtual_memory
from multiprocessing import cpu_count


# class Sample(object):
#     def __getattr__(self, name):
#         self.__dict__[name] = Sample()
#         return self.__dict__[name]


class Sample(object):
    def __int__(self, nanopore, illumina, assembly):
        self.nanopore = nanopore
        self.illumina = illumina
        self.assembly = assembly


class Nanopore(Sample):
    pass


class Illumina(Sample):
    pass


class Assembly(Sample):
    pass


class Bam(Sample):
    pass


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
    def check_dependencies(output_folder):
        """
        python=3.10.11 nextpolish=1.4.1 bwa=0.7.17 samtools=1.17 \
        porechop=0.2.4 filtlong=0.2.1 minimap2=2.26 flye=2.9.2 shasta=0.11.1 qualimap=2.2.2d bbmap=39.01 bandage=0.8.1 \
        fastp=0.22.0 ntedit=1.3.5 polypolish=0.5.0 pandas=1.5.3 seqtk=1.4 quast=5.2.0 medaka=1.8.0 mummer4=4.0.0rc1 \
        gnuplot=5.4.5 plotly=5.15.0
        """
        Methods.make_folder(output_folder)
        log_file = output_folder + '/log.txt'

        # Not being used right now because versions are captured in the requirements.txt file
        with open(log_file, 'w') as f:
            # Python
            # Python 3.10.11  # Outputs just this line
            stdout = subprocess.check_output(['python', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('Python v{}\n'.format(version))

            # Seqtk
            # Version: 1.4-r122  # Got to fetch this line
            try:
                stdout = subprocess.check_output(['seqtk'], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                stdout = e.output
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('Version:'):
                    version = line.split(' ')[1]
                    f.write('Seqtk v{}\n'.format(version))
                    break

            # Quast
            # QUAST v5.2.0  # Only output this line
            stdout = subprocess.check_output(['quast', '--version'], stderr=subprocess.STDOUT)
            f.write('{}\n'.format(stdout.decode('utf-8').rstrip()))

            # medaka
            # medaka 1.8.0  # only output this line
            stdout = subprocess.check_output(['medaka', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('Medaka v{}\n'.format(version))

            # mummer
            # 4.0.0rc1  # only output this line
            stdout = subprocess.check_output(['mummer', '--version'], stderr=subprocess.STDOUT)
            f.write('Mummer v{}\n'.format(stdout.decode('utf-8').rstrip()))

            # Nextpolish
            # nextPolish 1.4.1  # only this line output
            stdout = subprocess.check_output(['nextPolish', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('NextPolish v{}\n'.format(version))

            # Flye
            # 2.9.2-b1786  # only output this line
            stdout = subprocess.check_output(['flye', '--version'], stderr=subprocess.STDOUT)
            f.write('Flye v{}\n'.format(stdout.decode('utf-8').rstrip()))

            # Shasta
            # Shasta Release 0.11.1  # Oly output this line
            stdout = subprocess.check_output(['shasta', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[2].rstrip()
            f.write('Shasta v{}\n'.format(version))

            # Qualimap
            # QualiMap v.2.2.2-dev  # Got to fetch this line
            stdout = subprocess.check_output(['qualimap', '--version'], stderr=subprocess.STDOUT)
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('QualiMap v'):
                    version = line.split(' ')[1]
                    f.write('Qualimap v{}\n'.format(version))
                    break

            # Bandage
            # Version: 0.8.1  # only output this line
            stdout = subprocess.check_output(['Bandage', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('Bandage v{}\n'.format(version))

            # Fastp
            # fastp 0.22.0  # only output this line
            stdout = subprocess.check_output(['fastp', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('Fastp v{}\n'.format(version))

            # ntedit
            # ntedit version 1.3.5  # Got to fetch this line
            stdout = subprocess.check_output(['ntedit', '--version'], stderr=subprocess.STDOUT)
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('ntedit version'):
                    version = line.split(' ')[2]
                    f.write('ntEdit v{}\n'.format(version))
                    break

            # Polypolish
            # Polypolish v0.5.0  # only output this line
            stdout = subprocess.check_output(['polypolish', '--version'], stderr=subprocess.STDOUT)
            version = stdout.decode('utf-8').split(' ')[1].rstrip()
            f.write('Polypolish {}\n'.format(version))

            # Porechop
            # 0.2.4  # Outputs just this line
            stdout = subprocess.check_output(['porechop', '--version'], stderr=subprocess.STDOUT)
            f.write('Porechop v{}\n'.format(stdout.decode('utf-8').rstrip()))

            # BBtools suite
            # BBMap version 39.01  # Got to fetch that line from output
            stdout = subprocess.check_output(['bbduk.sh', '--version'], stderr=subprocess.STDOUT)
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('BBMap version'):
                    version = line.split(' ')[2]
                    f.write('BBtools v{}\n'.format(version))
                    break

            # Filtlong
            stdout = subprocess.check_output(['filtlong', '--version'], stderr=subprocess.STDOUT)
            # Filtlong v0.2.1  # Outputs just this line
            f.write('{}\n'.format(stdout.decode('utf-8').rstrip()))

            # SAMtools
            # samtools 1.17  # Got to fetch that line from output
            stdout = subprocess.check_output(['samtools', 'version'], stderr=subprocess.STDOUT)
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('samtools'):
                    version = line.split(' ')[1]
                    f.write('SAMtools v{}\n'.format(version))
                    break

            # BCFtools
            # BCFtools 1.17  # Got to fetch that line from output
            try:
                stdout = subprocess.check_output(['bcftools', 'version'], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                raise Exception("BCFtools not properly installed.")
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('bcftools'):
                    version = line.split(' ')[1]
                    f.write('BCFtools v{}\n'.format(version))
                    break

            # Minimap2
            # 2.26-r1175  # Outputs just this line
            stdout = subprocess.check_output(['minimap2', '--version'], stderr=subprocess.STDOUT)
            f.write('Minimap2 v{}\n'.format(stdout.decode('utf-8').rstrip()))

            # BWA
            # Version: 0.7.17-r1188  # Got to fetch that line from output
            try:
                stdout = subprocess.check_output(['bwa'], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                stdout = e.output
            for line in stdout.decode('utf-8').split('\n'):
                if line.startswith('Version:'):
                    version = line.split(' ')[1]
                    f.write('BWA v{}\n'.format(version))
                    break

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

                    # Initiate objects
                    sample_obj = Sample()
                    sample_obj.nanopore = Nanopore()
                    sample_obj.illumina = Illumina()  # Prep dictionary for Illumina, if present
                    sample_obj.assembly = Assembly()  # Prep dictionary for assembly
                    sample_obj.bam = Bam()  # Prep dictionary for mapping

                    # Add path to object
                    sample_dict[sample] = sample_obj
                    sample_dict[sample].nanopore.raw = file_path
                    sample_dict[sample].illumina.raw = list()
                    sample_dict[sample].illumina.trimmed = list()
                    sample_dict[sample].assembly.raw = ''
                    sample_dict[sample].bam.raw = ''

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

                    # Initiate objects
                    sample_obj = Sample()
                    sample_obj.illumina = Illumina()

                    if '_R1' in filename:
                        sample_dict[sample].illumina.raw.insert(0, file_path)
                    elif '_R2' in filename:
                        sample_dict[sample].illumina.raw.insert(1, file_path)

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
    def format_fasta(input_fasta, output_fasta):
        fasta_dict = dict()

        # Parse fasta file into dictionary
        with open(input_fasta, 'r') as in_handle:
            header = ''
            seq = list()  # Use list to store sequence information
            for line in in_handle:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seq:
                        fasta_dict[header] = ''.join(seq)  # Store in dictionary
                        seq = list()  # Empty seq
                    header = line
                else:
                    seq.append(line)
                # For the last entry
                fasta_dict[header] = ''.join(seq)

        # Write to file with a width of 80 character max for the sequence
        # Any additional sequences will be written of the next line
        with open(output_fasta, 'w') as out_handle:
            for title, seq in fasta_dict.items():
                out_handle.write('{}\n{}\n'.format(title, '\n'.join(textwrap.wrap(seq, 80, break_long_words=True))))


    @staticmethod
    def fix_mummerplot():
        # Fix gnuplot path in mummerplot
        mummerplot_path = shutil.which('mummerplot')
        gnuplot_path = shutil.which('gnuplot')

        # Check if file needs to be modified
        bad = False
        with open(mummerplot_path, 'r') as f:
            if '$GNUPLOT_EXE = \'false\'' in f.read():
                bad = True

        if bad:
            with open(mummerplot_path + '.tmp', 'w') as out_fh:
                with open(mummerplot_path, 'r') as in_fh:
                    for line in in_fh:
                        if '$GNUPLOT_EXE = \'false\'' in line:
                            line = 'my $GNUPLOT_EXE = \'{}\';'.format(gnuplot_path)
                            out_fh.write(line)
                        else:
                            out_fh.write(line)
            shutil.move(mummerplot_path + '.tmp', mummerplot_path)
            os.chmod(mummerplot_path, 0o0777)
