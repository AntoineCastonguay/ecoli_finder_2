import subprocess
import os
import sys
from concurrent import futures
import pathlib
from shutil import move
from psutil import virtual_memory
from multiprocessing import cpu_count
import gzip
from itertools import groupby
from glob import glob
import pysam
import warnings


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        if not os.path.exists(my_input):
            raise Exception('Please select an existing file or folder as input.')

        # Check if folder
        if os.path.isdir(my_input):
            file_list = os.listdir(my_input)  # List content of folder
        else:  # elif os.path.isfile(my_input):
            file_list = [my_input]

        # if folder is not empty and all files have the accepted extensions
        if not all([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_ref(ref):
        if not os.path.isfile(ref):
            raise Exception('The reference file provided does not exist')

        with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as f:
            first_header = f.readline()
            first_character = first_header[0]
            if first_character != '>':
                raise Exception('The reference file provided does not appear to be a valid fasta file.')

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
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    
    @staticmethod
    def get_sam_file(folder):
        for filename in os.listdir(folder):
            if filename.endswith('.sam'):
                return os.path.join(folder, filename)
        return None

    @staticmethod
    def list_files_in_folder(folder, extensions):
        if isinstance(extensions, str):
            extensions = [extensions]
        files = []
        for ext in extensions:
            files.extend(glob(f"{folder}/*.{ext}"))
        return files

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass
        
    @staticmethod
    def alignment(genome, primer, output):
        # Récupère les fichiers dans le dossier primer avec les extensions spécifiées
        extensions = ['fa', 'fasta']
        file = Methods.list_files_in_folder(primer, extensions)
        mon_dict = {}

        for f in file:
            name_file = os.path.basename(f)
            name = os.path.splitext(name_file)[0]
            if "_" in name:
                base, orientaition = name.rsplit('_', 1)
                orientaition = orientaition.upper()
                if orientaition != "R" and orientaition != "F":
                    print('Error : name of primer file incorrect need _F and _R. (primerX_F.fa and primerX_R.fa)')
            else:
                print('Error : name of primer file incorrect need _F and _R . (primerX_F.fa and primerX_R.fa)')
            if base not in mon_dict:
                mon_dict[base] = {}
            mon_dict[base][orientaition] = f

        # Crée le dossier de sortie si nécessaire
        os.makedirs(output, exist_ok=True)

        # Alignment BWA
        BWA_index_cmd = ['bwa', 'index', genome]
        subprocess.run(BWA_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Parcourt chaque entrée dans le dictionnaire pour effectuer les alignements
        for key, sub_dict in mon_dict.items():
            BWA_cmd = ['bwa', 'mem', genome, sub_dict['F'], sub_dict['R']]
            with open(f'{output}/BWA_output.sam', 'w') as outfile, open(os.devnull, 'w') as errfile:
                subprocess.run(BWA_cmd, stdout=outfile, stderr=errfile)

    @staticmethod
    def extract_primer_positions(sam_file, essentiel_file):
        primer_positions = {}
        with open(essentiel_file,'r') as essentiel:
            liste_gene_essentiel = [ligne.strip() for ligne in essentiel]

        with open(sam_file, 'r') as file:
            for line in file:
                # Ignorer les lignes d'en-tête
                if line.startswith('@'):
                    continue

                # Séparer la ligne en colonnes
                columns = line.strip().split('\t')
                read_id = columns[0]  # ID de la lecture
                flag = int(columns[1])  # Flag de la lecture
                position = int(columns[3])  # Position d'alignement
                qualite = columns[5]  #qualité alignment
                postion_mate = int(columns[7])
                length = int(columns[8])
                seq = columns[9]
                
                parties = read_id.split('-')
                id = "-".join(parties[:2])
                gene = parties[2]
                init = parties[3]

                if read_id in liste_gene_essentiel:
                    essentiel_gene = True
                else:
                    essentiel_gene = False

                list_var = [gene,position,postion_mate,length,qualite,essentiel_gene,init]

                if id not in primer_positions:
                    primer_positions[id] = {}
                primer_positions[id][flag] = list_var

        return primer_positions
    
    @staticmethod
    def write_result(data,output):
        print('Creation of result file...')

        Methods.make_folder(output)
        with open(f'{output}/output.txt', 'w') as f:
            f.write(f"id\tgene\tflag\tfirst_pos\tsecond_pos\tlength\tquality\tessentiel\tinit\n")
            for read_id, sub_dict in data.items():
                for flag, list_var in sub_dict.items():
                    f.write(f"{read_id}\t{list_var[0]}\t{flag}\t{list_var[1]}\t{list_var[2]}\t{list_var[3]}\t{list_var[4]}\t{list_var[5]}\t{list_var[6]}\n")