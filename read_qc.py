import os
import subprocess
from concurrent import futures
from bga_methods import Methods


class ReadQC(object):
    @staticmethod
    def read_qc_falco(sample, path_list, output_folder, cpu):
        for read_file in path_list:
            cmd = ['falco',
                   '--outdir', output_folder,
                   '--threads', str(cpu),
                   '-report-filename', output_folder + os.path.basename(read_file).split('.')[0] + '.html',
                   '-skip-data',
                   '-skip-summary',
                   read_file]

            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def read_qc_falco_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path_list, output_folder, int(cpu / parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: ReadQC.read_qc_falco(*x), args):
                pass

    @staticmethod
    def run_pycoqc(basecalled_folder, report_folder):
        # Need the "sequencing_summary.txt" file
        Methods.make_folder(report_folder)
        cmd = ['pycoQC',
               '-f', basecalled_folder + 'sequencing_summary.txt',
               '-o', report_folder + 'pycoQC_output.html']
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
