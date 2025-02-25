import os
from typing import OrderedDict
from Bio import Entrez
import pandas as pd
import gzip
import argparse 
import pathlib
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess as sb
import glob
import shutil
import re
import logging
import csv
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
import numpy as np
from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw
from Bio.Align import PairwiseAligner
import tempfile
import re
try:
    import pkg_resources
except ImportError:
    pkg_resources = None

def validate_email(value):
    if not re.match(r"[^@]+@[^@]+\.[^@]+", value):
        raise argparse.ArgumentTypeError(f"Invalid email address: {value}")
    return value

parser = argparse.ArgumentParser(description="This command-line tool facilitates the bulk downloading of various genomes, performs extraction of 16S and hypervariable regions, and processes the files to be suitable as inputs for the Qiime2 Klebotoxin library.")

parser.add_argument("-o", "--output", type=pathlib.Path,
                    help="Specifies the path for output files. If not provided, the tool will generate the necessary files in the current directory, creating a separate download_file for each downloaded genome and an output file for the produced files.")

parser.add_argument("-m","--list_genomes", nargs="+", type=str, 
                    default=["Klebsiella oxytoca", "Klebsiella variicola", "Klebsiella michiganensis","Klebsiella grimontii","Klebsiella aerogenes","Klebsiella pneumoniae"], 
                    help="Defines a list of genomes to process.")

parser.add_argument("--refseq", action="store_true", help="Download only RefSeq assemblies if set.")

parser.add_argument("-g","--genome_number", default="20", type=int,
                    help="Specifies the number of genomes to download.")

parser.add_argument("-d","--download_skip",action="store_true",
                    help="If specified, the tool will skip the genome downloading process. Useful if you already have the genomes and need to modify other parameters.")

parser.add_argument("-p","--primers_off", action="store_false",
                    help="If specified, the tool uses hardcoded nucleotide locations for the search of the hypervariable regions.")

parser.add_argument("-2","--pairwise2",action="store_true",
                    help="If specified, the tool uses the outdated pairwise2 approach. This will trigger a warning as it's not recommended due to deprecation. However, if selected, the tool will still provide the results.")

parser.add_argument("-s","--extraction_skip",action="store_true",
                    help="Skips the extraction of 16S and hypervariable regions. Useful if the only requirement is to convert files to '.qza' format.")

parser.add_argument("-c","--convert_qza",action="store_true",
                    help="If specified, the tool converts the output into Qiime artifacts (.qza) which can be used as input for klebotoxin.")

parser.add_argument("-e", "--email", help="Email address is required for NCBI's Entrez databases when accessed via Biopython", type=validate_email)

parser.add_argument("-a", "--api_key", help="API key for NCBI's Entrez databases when accessed via Biopython. Code can run without one but highly suggested to pass one especially if downloading large amount of genomes",
                     type=str)

parser.add_argument("-b","--blaoxy_check",action="store_true",
                    help="If specified, the tool does a blaoxy check verification.")

parser.add_argument("-k", "--consensus", action="store_true", help="Create a consensus sequence from extracted 16S rRNA sequences.")

parser.add_argument("-t","--standalone", action="store_true", help="Run in standalone mode without KoSCAPEbio plugin")

parser.add_argument("-x", "--bla_oxy", type=pathlib.Path,
                        help="Specifies the path for bla_oxy database. Only use in standalone mode if bla_database is located elsewhere")

parser.add_argument("-l", "--logging", action="store_true", help="Enable logging")

args = parser.parse_args()





def main():
    if args.logging:
        def setup_loggers(log_path, info_log_name, error_log_name):       
            os.makedirs(log_path, exist_ok=True)
            info_logger = logging.getLogger(info_log_name)
            info_file_handler = logging.FileHandler(os.path.join(log_path, f"{info_log_name}.log"))
            info_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            info_file_handler.setFormatter(info_formatter)
            info_logger.addHandler(info_file_handler)
            info_logger.setLevel(logging.INFO)
            error_logger = logging.getLogger(error_log_name)
            error_file_handler = logging.FileHandler(os.path.join(log_path, f"{error_log_name}.log"))
            error_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            error_file_handler.setFormatter(error_formatter)
            error_logger.addHandler(error_file_handler)
            error_logger.setLevel(logging.ERROR)

            return info_logger, error_logger
                
    if not args.download_skip and args.email is None:
        parser.error("--email is needed when downloading genomes from NCBI's Entrez databases using Biopython. Please provide a valid email address.")
    
    if args.bla_oxy and not args.standalone:
        parser.error("--bla_oxy can only be used with --standalone")

    if args.output==None:
        print("Missing the output path, defaulting to location of the script")
        output=str(pathlib.Path.cwd())
    else:
        output=str(args.output)

    if not os.path.exists(f"{output}/Genomes"):
        os.makedirs(f"{output}/Genomes")
    if not os.path.exists(f"{output}/Output_files"):
        os.makedirs(f"{output}/Output_files")
    log_directory = os.path.join(f"{output}/PipelineResultsArchive/Logs")
    csv_directory = os.path.join(f"{output}/PipelineResultsArchive/CSV")
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    if not os.path.exists(csv_directory):
        os.makedirs(csv_directory)
    main_info_path=f"{output}/PipelineResultsArchive"
    log_path=f"{output}/PipelineResultsArchive/Logs"
    csv_path=f"{output}/PipelineResultsArchive/CSV"


    genome_path=f"{output}/Genomes"
    output_path=f"{output}/Output_files"
    

    organisms=args.list_genomes
    max_genomes_per_organism = int(args.genome_number)

    if max_genomes_per_organism < 10:
        retmax_ = max_genomes_per_organism + 20
    elif 10 <= max_genomes_per_organism < 50:
        retmax_ = max_genomes_per_organism + 50
    elif 50 <= max_genomes_per_organism < 100:
        retmax_ = max_genomes_per_organism + 100
    else:
        retmax_ = max_genomes_per_organism + max_genomes_per_organism
    if not args.download_skip:
        fasta_path = f"{genome_path}/fasta"
        blast_results_path = f"{genome_path}/blast_results"

        if not os.path.exists(fasta_path):
            os.makedirs(fasta_path)

        if not os.path.exists(blast_results_path):
            os.makedirs(blast_results_path)        
        def is_bla_oxy_present(fasta_path, blast_db):
            blast_output_path = f"{blast_results_path}/{os.path.basename(fasta_path)}.out"
            print(f"Running BLAST on {fasta_path} against {blast_db}...")
            
            blast_command = [
                'blastn',
                '-query', fasta_path,
                '-db', blast_db,
                '-evalue', '0.04',
                '-outfmt', '6',
                '-out', blast_output_path
            ]

            try:
                sb.run(blast_command, check=True, capture_output=True)
                print("BLAST search completed successfully.")
            except sb.CalledProcessError as e:
                print(f"BLAST search failed with error: {e.stderr.decode()}")
                return None, None

            with open(blast_output_path, 'r') as f:
                line = f.readline()
                if line:
                    fields = line.strip().split('\t')
                    top_hit_species = fields[1]
                    percent_identity = fields[2]
                    top_hit_species = re.sub(r'_\d+$', '', top_hit_species)
                    
                    print(f"Top BLAST hit species: {top_hit_species.replace('_', ' ')}")
                    print(f"Percent identity: {percent_identity}%")
                    
                    return top_hit_species.replace("_", " "), percent_identity
                else:
                    print(f"No BLAST hits found for {fasta_path}")
            
            return None, None
        
        
        Entrez.email = args.email 
        if args.api_key==None:
            api_key=""
            print("Warning missing argument '--api_key' the download may be slow or even stop. High suggestion to pass the API key")
        else:
            api_key=str(args.api_key)

        Entrez.api_key = api_key
        def download_genbank_assembly(ftp_path, genome_path):
            assembly_name = ftp_path.split("/")[-1]
            filename = f"{assembly_name}_genomic.gbff.gz"
            full_ftp_path = f"{ftp_path}/{filename}"
            full_https_path = full_ftp_path.replace("ftp://", "https://") 
            print(f"Downloading: {filename}")   
            os.system(f"wget {full_https_path} -c -t 10 -P {genome_path}")
            return assembly_name  
        def get_species(description):
            species = description.split()[0]

            step_1=species.replace('_', ' ')
            step_2=step_1.replace(' 1','')
            step_3=step_2.replace(' 5','')
            return step_3

        def verify_and_correct_species_name(gbff_file_path, expected_species_name, output_path):
            change_log_file = os.path.join(output_path, "renamed_genomes.csv")

            with gzip.open(gbff_file_path, 'rt') as file:
                content = file.readlines()

            changes_made = {}
            updated_content = []
            occurrences = 0

            for line in content:
                original_line = line
                if line.startswith("SOURCE") or line.startswith("DEFINITION") or line.startswith("LOCUS"):
                    current_species_name = re.search(r'(\b[A-Z][a-z]+ [a-z]+\b)', line)
                    if current_species_name:
                        current_species_name = current_species_name.group()
                        if current_species_name.lower() != expected_species_name.lower():
                            line = re.sub(re.escape(current_species_name), expected_species_name, line, flags=re.IGNORECASE)
                            changes_made[original_line.strip()] = line.strip()
                            occurrences += 1
                updated_content.append(line)

            if occurrences > 0:
                print(f"Checked genome {os.path.basename(gbff_file_path)}: Naming convention failed. Correcting names. \n")

                uncompressed_path = gbff_file_path.replace('.gz', '')
                with open(uncompressed_path, 'w') as file:
                    file.writelines(updated_content)

                with open(uncompressed_path, 'rb') as f_in:
                    with gzip.open(gbff_file_path, 'wb') as f_out:
                        f_out.writelines(f_in)

                os.remove(uncompressed_path)
                if not os.path.exists(change_log_file):
                    with open(change_log_file, 'w', newline='') as file:
                        writer = csv.writer(file)
                        writer.writerow(["file_name", "old_definition", "old_source", "old_locus", "corrected_name", "number_of_corrections"])

                with open(change_log_file, 'a', newline='') as file:
                    writer = csv.writer(file)
                    for old_line, new_line in changes_made.items():
                        writer.writerow([os.path.basename(gbff_file_path), old_line, "", "", new_line, occurrences])

            else:
                print(f"Checked genome {os.path.basename(gbff_file_path)}: Naming convention passed.")
            
            return gbff_file_path



        if args.standalone:
            if args.bla_oxy:
                blaoxy_fasta_path = args.bla_oxy / "blaoxy_database.fasta" #weird 
                bla_db_path = args.bla_oxy / "blaoxy_db"
            else:
                blaoxy_fasta_path = pathlib.Path("bla_database/blaoxy_database.fasta")
                bla_db_path = pathlib.Path("bla_database/blaoxy_db")
            blaoxy_organisms = [get_species(record.description) for record in SeqIO.parse(blaoxy_fasta_path, "fasta")]
        else:
            blaoxy_fasta_path = pkg_resources.resource_filename('q2_koscapebio.database_curation', 'bla_database/blaoxy_database.fasta')
            blaoxy_organisms = [get_species(record.description) for record in SeqIO.parse(blaoxy_fasta_path, "fasta")]


        csv_file_path = os.path.join(csv_path, 'genome_aquiring.csv')

        with open(csv_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Specie Name', 'File Name', 'File Definition', 'BlaOxy Check','Percent identitiy',  'Saved?'])
        for organism in organisms:
            from io import BytesIO
            print(f"Working on {organism}")
            assembly_data = []
            handle = Entrez.esearch(db="assembly", term=organism + "[organism]", retmax=str(retmax_))
            data = handle.read() 
            record = Entrez.read(BytesIO(data), validate=False) 
            assembly_ids = record["IdList"]
            def get_assembly_summary(assembly_id):
                try:
                    handle = Entrez.esummary(db="assembly", id=assembly_id)
                    response_content = handle.read()
                    handle.close()
                    record = Entrez.read(BytesIO(response_content), validate=False)
                    return record
                except Exception as e:
                    print(f"An error occurred while fetching assembly summary: {e}")
                    print(f"The response content was: {response_content}")
                    raise
            for assembly_id in assembly_ids:
                try:
                    record = get_assembly_summary(assembly_id)
                    if not record["DocumentSummarySet"]["DocumentSummary"]:
                        print(f"No summary found for assembly ID {assembly_id}.")
                        continue
                    if args.refseq:
                        if record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'].startswith('GCF'):
                            assembly_accession = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
                            ftp_path = record["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_RefSeq")
                    else:
                        assembly_accession = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
                        ftp_path = record["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_GenBank")
                except IndexError as e:
                    print(f"An error occurred: {e}")
                    print(f"The record retrieved was: {record}")
                    continue
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    sys.exit(1)
                
                if not ftp_path:
                    print(f"No FTP path found for assembly id {assembly_id}. Skipping this assembly.")
                    continue
                
                assembly_data.append((assembly_accession, ftp_path))
            
            count = 0
            if args.logging:
                info_logger, error_logger = setup_loggers(log_path, 'conversion_info', 'conversion_errors')

            for assembly_accession, ftp_path in assembly_data:
                assembly_name = ftp_path.split("/")[-1]
                filename = f"{assembly_name}_genomic.gbff.gz"
                gzipped_gbff = os.path.join(genome_path, filename)
                if os.path.exists(gzipped_gbff):
                    print(f"{filename} already exists. Skipping download.")
                    if args.logging:
                        info_logger.info(f"{filename} already exists. Skipping download.")
                    continue
                else:
                    download_genbank_assembly(ftp_path, genome_path)
                organism_formatted_for_comparison = organism.replace("_", " ")
                renamed_gbff_path = verify_and_correct_species_name(gzipped_gbff, organism_formatted_for_comparison, output_path)
                try:
                    fasta_file = gzipped_gbff.replace(".gbff.gz", ".fasta")
                    SeqIO.convert(gzip.open(gzipped_gbff, 'rt'), "genbank", fasta_file, "fasta")
                    if args.logging:
                        info_logger.info(f"Successfully converted {gzipped_gbff} to {fasta_file}.")

                except Exception as e:
                    if args.logging:
                        error_logger.error(f"Failed to convert {gzipped_gbff} to {fasta_file} due to: {str(e)}")
                        error_logger.info(f"Removed corrupted file {gzipped_gbff}.")
                    os.remove(gzipped_gbff)
                    continue 



                should_do_blast_check = organism_formatted_for_comparison in blaoxy_organisms and args.blaoxy_check

                if should_do_blast_check:
                    print(f"Performing bla-oxy check for {organism}...")
                    if args.standalone:
                        percent_identity=0
                        blast_top_hit_species, percent_identity = is_bla_oxy_present(fasta_file, bla_db_path)
                    else:
                        percent_identity=0
                        bla_db_path = pkg_resources.resource_filename('q2_koscapebio.database_curation', 'bla_database/blaoxy_db')
                        blast_top_hit_species, percent_identity = is_bla_oxy_present(fasta_file, bla_db_path)

        
            

                    if blast_top_hit_species != organism.replace("_"," "):
                        blast_results=f"Failed - Top result was {blast_top_hit_species}"
                        file_saved="No"
                        print(f"BLA-OXY check failed for {organism}. Removing files and renaming BLAST output. \n")
                        blast_output_file = f"{blast_results_path}/{os.path.basename(fasta_file)}.out"
                        new_blast_output_file = f"{blast_results_path}/{os.path.basename(fasta_file)}_REMOVED_.out"
                        os.rename(blast_output_file, new_blast_output_file)

                        os.remove(gzipped_gbff)
                        os.remove(fasta_file)
                        continue
                    else:
                        print(f"BLA-OXY check passed for {organism}. Moving fasta file to {fasta_path}/{os.path.basename(fasta_file)} \n")
                        blast_results=f"Passed - Top result was {blast_top_hit_species}"
                        file_saved="Yes"
                        shutil.move(fasta_file, f"{fasta_path}/{os.path.basename(fasta_file)}")
                        count += 1
                else:
                    percent_identity=0
                    print(f"No bla-oxy check required for {organism}. Moving fasta file to {fasta_path}/{os.path.basename(fasta_file)} \n")
                    blast_results=f"No bla-oxy check required for {organism_formatted_for_comparison}"
                    file_saved="Yes"
                    shutil.move(fasta_file, f"{fasta_path}/{os.path.basename(fasta_file)}")
                    count += 1



                with gzip.open(renamed_gbff_path, "rt") as handle:
                    for record in SeqIO.parse(handle, "genbank"):
                        definition = record.description 
                        break     
           
                with open(csv_file_path, mode='a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow([organism_formatted_for_comparison, filename, definition, blast_results,percent_identity, file_saved])

                if count >= max_genomes_per_organism:
                    print(f"Reached maximum number of genomes for {organism}. Breaking loop. \n")
                    break

    def remove_file_if_exists(filepath):
        if os.path.isfile(filepath):
            os.remove(filepath)

    output_files_main = ["alignments.txt"]
    output_files_fasta = ["16S.fasta", "V3_V4_clean.fasta", 
                    "V3_V4_duplicated.fasta", "V4_V5_clean.fasta", "V4_V5_duplicated.fasta"]
    output_files_qza = ["16S.qza", "V3_V4_clean.qza", 
                    "V3_V4_duplicated.qza", "V4_V5_clean.qza", "V4_V5_duplicated.qza"]

    directories_files = {"": output_files_main, "FASTA": output_files_fasta, "QZA": output_files_qza}

    for dir, files in directories_files.items():
        for file in files:
            remove_file_if_exists(os.path.join(output_path, dir, file))


    def format_fasta_header(record,counter):
        counter=counter+1
        try:
            description_parts = record.description.split(" ")
            if 'strain' in description_parts:
                strain_index = description_parts.index('strain')
                species_strain_name = "_".join(description_parts[:strain_index + 2])
            else:
                species_strain_name = "_".join(description_parts[:4])
            species_strain_name = species_strain_name.replace(",", "").replace("(", "").replace(")", "")
            accession_id = record.id.split(".")[0]
            new_header = f"{species_strain_name}_{accession_id}_gene_copy-{counter}"
        except Exception as e:
            print(f"Error formatting header for record {record.id}: {e}")
            new_header = "Unknown_header"
        return new_header

    def align_sequences(sequences):
        """Align sequences using MAFFT and return the alignment object, with detailed logging."""
        try:
            if args.logging:
                info_logger, error_logger = setup_loggers(log_path, "mafft_info", "mafft_errors")
            with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta") as input_fasta:
                SeqIO.write(sequences, input_fasta, 'fasta')
                input_fasta.flush()
                with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta") as output_fasta:
                    process = sb.run(['mafft', '--auto', input_fasta.name], stdout=output_fasta, stderr=sb.PIPE, text=True)
                    output_fasta.seek(0)
                    alignment = AlignIO.read(output_fasta, 'fasta')
                    if args.logging:
                        info_logger.info("MAFFT alignment completed successfully.")
                    
                    if process.stderr:
                        print("MAFFT reported the following error:\n" + process.stderr)
                        if args.logging:
                            error_logger.error("MAFFT reported the following error:\n" + process.stderr)
                        
            return alignment
        except Exception as e:
            error_logger.error("MAFFT alignment failed: " + str(e))
            return None

    def create_consensus_sequence_new(alignment, consensus_threshold):
        """Generate consensus sequences from groups within an alignment based on similarity and assign sequential consensus IDs."""
        similarity_threshold=99.9
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        
        groups = []
        processed_sequences = set()

        for i in range(len(alignment)):
            if i in processed_sequences:
                continue
            current_group = [alignment[i]]
            for j in range(i + 1, len(alignment)):
                if j in processed_sequences:
                    continue

                seq1 = alignment[i].seq
                seq2 = alignment[j].seq
                score = aligner.score(seq1, seq2)
                max_possible_score = min(len(seq1), len(seq2)) * aligner.match_score
                similarity = (score / max_possible_score) * 100

                if similarity >= similarity_threshold:
                    current_group.append(alignment[j])
                    processed_sequences.add(j)

            groups.append(current_group)
            processed_sequences.add(i)

        consensus_sequences = []
        group_index=1
        for group in groups:
            consensus_seq = ""

            for i in range(len(group[0].seq)):  
                base_count = {}
                for record in group:
                    base = record.seq[i]
                    base_count[base] = base_count.get(base, 0) + 1

                total = sum(base_count.values())
                consensus_base = None
                for base, count in base_count.items():
                    if count / total >= consensus_threshold:
                        consensus_base = base
                        break

                if consensus_base:
                    consensus_seq += consensus_base
                else:
                    consensus_seq += 'N'  

            consensus_id = re.sub(r'_gene_copy-\d+', '', group[0].id) 
            consensus_id = f"{consensus_id}_consensus-{group_index}"
            group_index += 1  

            consensus_record = SeqRecord(Seq(consensus_seq), id=consensus_id, description="")
            consensus_sequences.append(consensus_record)

        return consensus_sequences, groups

    def clean_sequence(seq_record):
        cleaned_seq = str(seq_record.seq).upper()
        cleaned_seq_record = SeqRecord(Seq(cleaned_seq), id=seq_record.id, description=seq_record.description)
        return cleaned_seq_record

    def perform_similarity_check(temp_fasta_path, similarity_threshold=99.9):
        seq_records = list(SeqIO.parse(temp_fasta_path, "fasta"))
        if len(seq_records) < 2:
            return seq_records, [] 
        
        keep = []
        outliers = []
        processed_sequences = set() 
        aligner = PairwiseAligner()
        aligner.mode = 'local' 

        for i in range(len(seq_records)):
            if i in processed_sequences:
                continue
            
            found_similar = False
            group_to_process = [i]
            
            for j in range(i + 1, len(seq_records)):
                if j in processed_sequences:
                    continue

                seq1 = seq_records[i].seq
                seq2 = seq_records[j].seq
                score = aligner.score(seq1, seq2)
                max_possible_score = min(len(seq1), len(seq2)) * aligner.match_score
                similarity = (score / max_possible_score) * 100
                if similarity >= similarity_threshold:
                    found_similar = True
                    group_to_process.append(j)
                    processed_sequences.add(j)
            if found_similar:
                for seq_idx in group_to_process:
                    keep.append(seq_records[seq_idx])
                    processed_sequences.add(seq_idx)
            else:
                outliers.append(seq_records[i])
        
        return keep, outliers


        
    def extract_16S_rRNA_sequences(genome_path, temp_fasta_path, rRNA_regex, csv_writer, error_logger):
        with gzip.open(genome_path, "rt") as handle, open(temp_fasta_path, "w") as temp_out_handle:
            counter = 0
            rRNA_count = 0
            short_sequences = 0
            definition_rna = ""

            try:
                for seq_record in SeqIO.parse(handle, "genbank"):
                    if not definition_rna:
                        definition_rna = seq_record.description

                    for feature in seq_record.features:
                        if feature.type == "rRNA":
                            rRNA_count += 1
                            product = feature.qualifiers.get("product", [""])[0]
                            if rRNA_regex.search(product):
                                if len(feature) >= 1400:
                                    extracted_seq = feature.extract(seq_record)
                                    new_record = SeqIO.SeqRecord(
                                        extracted_seq.seq,
                                        id=format_fasta_header(seq_record, counter),
                                        description=""
                                    )
                                    SeqIO.write(new_record, temp_out_handle, "fasta")
                                    counter += 1
                                else:
                                    short_sequences += 1

                return counter, short_sequences, rRNA_count, definition_rna

            except Exception as e:
                error_logger.error(f"Failed to process file {genome_path}: {e}")
                return 0, 0, 0, ""

    if not args.extraction_skip:
        csv_file_path = os.path.join(csv_path, '16S_extraction.csv')
        with open(csv_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['File Name',"File Definition", 'Number of 16S Extracted', 'Reason for No Extraction', 'Number of Keep Sequences', 'Number of Outlier Sequences','Number of different sequences saved'])

        rRNA_regex = re.compile(r"16S\s*r(ibosomal)?\s*R(NA)?", re.IGNORECASE)
        if args.logging:
            info_logger, error_logger = setup_loggers(log_path, '16S-extraction_info', '16S-extraction_errors')

        with open(csv_file_path, mode='a', newline='') as file:
            csv_writer = csv.writer(file)

            for file in os.listdir(genome_path):
                if file.endswith(".gbff.gz"):
                    print(f"Working on the file {file}\n")
                    temp_fasta_file = tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix=".fasta")
                    renamed_gbff_path = os.path.join(genome_path, file)
                    
                    num_extracted, short_sequences, rRNA_count, definition_rna = extract_16S_rRNA_sequences(renamed_gbff_path, temp_fasta_file.name, rRNA_regex, csv_writer, error_logger)

                    if num_extracted == 0:
                        reason = (
                            "No 16S rRNA sequences found" if rRNA_count == 0 else
                            "All rRNA sequences shorter than 1400 bases" if short_sequences > 0 else
                            "No matching 16S rRNA product criteria"
                        )
                        csv_writer.writerow([file,definition_rna, num_extracted, reason, 0, 0])
                    else:
                        if args.consensus:
                            keep, outliers = perform_similarity_check(temp_fasta_file.name)
                            groups=0
                            if not keep:
                                print("No sequences to align for consensus.")
                                csv_writer.writerow([file,definition_rna, num_extracted, 'No sequences for consensus alignment', 0, len(outliers),0])
                                continue
                            else:
                                if len(keep) == 1:
                                    cleaned_consensus_sequences = clean_sequence(keep[0])
                                    groups=1
                                elif len(keep) == 2:
                                    cleaned_consensus_sequences = []
                                    
                                    for idx, seq_record in enumerate(keep, start=1):
                                        new_id = re.sub(r'_gene_copy-\d+', '', seq_record.id)
                                        new_id = f"{new_id}_consensus-{idx}"
                                        
                                        renamed_seq_record = SeqRecord(seq_record.seq, id=new_id, description="")
                                        
                                        cleaned_consensus_sequence = clean_sequence(renamed_seq_record)
                                        cleaned_consensus_sequences.append(cleaned_consensus_sequence)

                                    groups = len(cleaned_consensus_sequences)
                                else:
                                    alignment = align_sequences(keep)
                                    if alignment:
                                        if len(keep) == 3:
                                            consensus_threshold = 0.6
                                        else:
                                            consensus_threshold = 0.7
                                        consensus_sequence , groups_list = create_consensus_sequence_new(alignment, consensus_threshold)
                                        groups=len(groups_list)
                                    cleaned_consensus_sequences = []
                                    for consensus_seq in consensus_sequence:
                                        cleaned_consensus_sequence = clean_sequence(consensus_seq)
                                        cleaned_consensus_sequences.append(cleaned_consensus_sequence)


                                with open(os.path.join(output_path, "removed-seq.fasta"), "a") as removed_output_handle:
                                    SeqIO.write(outliers, removed_output_handle, "fasta")
                                with open(os.path.join(output_path, "16S.fasta"), "a") as output_handle:
                                    SeqIO.write(cleaned_consensus_sequences, output_handle, "fasta")

                                csv_writer.writerow([file,definition_rna, num_extracted, '', len(keep), len(outliers),groups])

                        else:
                            with open(temp_fasta_file.name, "r") as temp_in_handle:
                                temp_contents = temp_in_handle.read()
                            with open(os.path.join(output_path, "16S.fasta"), "a") as final_output_handle:
                                final_output_handle.write(temp_contents)

                            csv_writer.writerow([file,definition_rna, num_extracted, '', 0, 0,num_extracted])
                    os.remove(temp_fasta_file.name)

        records = list(SeqIO.parse(f"{output_path}/16S.fasta", "fasta"))
        if args.primers_off==False:
            v3_v4_records = []
            v4_v5_records = []
            v3_v5_records = []
            v4_records = []
            for record in records:
                v3_v4_records.append(record[433:805])  # Assuming the V3 region starts at 433 and V4 ends at 805
                v4_v5_records.append(record[515:926])  # Assuming the V4 region starts at 515 and V5 ends at 926
                v4_records.append(record[515:805])  # Assuming the V4 region starts at 515 and ends at 805
                v3_v5_records.append(record[433:926])  # Assuming the V3 region starts at 433 and V5 ends at 926
        if args.primers_off==True:
            if args.pairwise2==True:
                from Bio import pairwise2
                from Bio.pairwise2 import format_alignment
                def find_primer_positions(primer, sequence,alignment_file):
                    alignments = pairwise2.align.localms(sequence, primer, 2, -1, -0.5, -0.1)
                    top_alignment = alignments[0]
                    
                    with open(alignment_file, "a") as file:
                        file.write(format_alignment(*top_alignment))
                    start = top_alignment[3]
                    end = top_alignment[4]
                    print(f"Alignment found for primer {primer} in sequence {sequence[:10]}... Start: {start}, End: {end}")
                    return start, end
            if args.pairwise2==False:
                def find_primer_positions(primer, sequence, alignment_file):
                    aligner = PairwiseAligner()
                    aligner.mode = 'local'
                    aligner.match_score = 2
                    aligner.mismatch_score = -1
                    aligner.open_gap_score = -0.5
                    aligner.extend_gap_score = -0.1

                    alignments = aligner.align(sequence, primer)
                    if not alignments:
                        print(f"No alignments found for primer {primer} in sequence {sequence[:10]}...")
                        return None, None

                    best_alignment = max(alignments, key=lambda x: x.score)
                    start, end = best_alignment.aligned[0][0]

                    with open(alignment_file, "a") as file:
                        file.write(str(best_alignment))

                    return start, end

            primer_341F = "CCTACGGGNGGCWGCAG"  
            primer_805R = "ATTAGAWACCCBNGTAGTCC"  

            primer_515F = "GTGCCAGCMGCCGCGGTAA"  
            primer_926R = "AAACTYAAAKGAATTGACGG"  

            primer_806R = "ATTAGAWACCCBDGTAGTCC"  
                           

            alignment_file =f"{output_path}/alignments.txt"
            v3_v4_records = []
            v4_v5_records = []    
            v4_records = []
            v3_v5_records=[]
            for record in records:
                seq_str = str(record.seq)
                
                start_341F, end_341F = find_primer_positions(primer_341F, seq_str, alignment_file)
                start_805R, end_805R = find_primer_positions(primer_805R, seq_str, alignment_file)
                start_515F, end_515F = find_primer_positions(primer_515F, seq_str, alignment_file)
                start_926R, end_926R = find_primer_positions(primer_926R, seq_str, alignment_file)
                start_806R, end_806R = find_primer_positions(primer_806R, seq_str, alignment_file)


                v3_v4_seq = record[end_341F:end_805R]
                v4_v5_seq = record[end_515F:end_926R]
                v4_seq = record[end_515F:start_806R]
                v3_v5_seq = record[end_341F:end_926R]

                if len(v3_v4_seq) <= 500:
                    v3_v4_records.append(v3_v4_seq)
                else:
                    print(f"V3-V4 sequence for record {record.id} was longer than 600bp and has been skipped")

                if len(v4_v5_seq) <= 500:
                    v4_v5_records.append(v4_v5_seq)
                else:
                    print(f"V4-V5 sequence for record {record.id} was longer than 600bp and has been skipped")
                
                if len(v4_seq) <= 350:
                    v4_records.append(v4_seq)
                else:
                    print(f"V4 sequence for record {record.id} was longer than 600bp and has been skipped")

                if len(v3_v5_seq) <= 800:
                    v3_v5_records.append(v3_v5_seq)
                else:
                    print(f"V3-V5 sequence for record {record.id} was longer than 600bp and has been skipped")


        SeqIO.write(v3_v4_records, f"{output_path}/V3_V4_duplicated.fasta", "fasta")
        SeqIO.write(v4_v5_records, f"{output_path}/V4_V5_duplicated.fasta", "fasta")
        SeqIO.write(v4_records, f"{output_path}/V4_duplicated.fasta", "fasta")
        SeqIO.write(v3_v5_records, f"{output_path}/V3_V5_duplicated.fasta", "fasta")

        shutil.copy(f"{output_path}/V3_V4_duplicated.fasta", f"{output_path}/V3_V4_clean.fasta")
        shutil.copy(f"{output_path}/V4_V5_duplicated.fasta", f"{output_path}/V4_V5_clean.fasta")
        shutil.copy(f"{output_path}/V4_duplicated.fasta", f"{output_path}/V4_clean.fasta")
        shutil.copy(f"{output_path}/V3_V5_duplicated.fasta", f"{output_path}/V3_V5_clean.fasta")

        def remove_duplicates(file):
            records = list(SeqIO.parse(file, "fasta"))
            records_df = pd.DataFrame([(record.id, record.description, str(record.seq)) for record in records], columns=["id", "description", "seq"])
            records_df.drop_duplicates(subset="seq", inplace=True) 
            deduplicated_records = [SeqRecord(Seq(seq), id=id, description=desc) for id, desc, seq in records_df.to_records(index=False)]
            SeqIO.write(deduplicated_records, file, "fasta")

        remove_duplicates(f"{output_path}/V3_V4_clean.fasta")
        remove_duplicates(f"{output_path}/V4_V5_clean.fasta")
        remove_duplicates(f"{output_path}/V4_clean.fasta")
        remove_duplicates(f"{output_path}/V3_V5_clean.fasta")
    if args.convert_qza==True:
        def run_qiime_import(input_path, output_path):
            cmd = [
                "qiime", "tools", "import", 
                "--input-path", input_path, 
                "--output-path", output_path, 
                "--type", "FeatureData[Sequence]"
            ]
            
            process = sb.Popen(cmd, stdout=sb.PIPE, stderr=sb.PIPE)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                print(f"An error occurred: {stderr.decode()}")
            else:
                print(f"Output: {stdout.decode()}")

        def convert_all_fasta_files(output_path):

            fasta_files = glob.glob(f"{output_path}/*_clean.fasta")

            for fasta_file in fasta_files:
                qza_file = os.path.splitext(fasta_file)[0] + ".qza"
                run_qiime_import(fasta_file, qza_file)

        convert_all_fasta_files(output_path)

    def organizing_files():
        fasta_dir = os.path.join(f"{output_path}/", "FASTA")
        qza_dir = os.path.join(f"{output_path}/", "QZA")

        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(qza_dir, exist_ok=True)

        for fasta_file in glob.glob(f"{output_path}/*.fasta"):
            shutil.move(fasta_file, fasta_dir)

        for qza_file in glob.glob(f"{output_path}/*.qza"):
            shutil.move(qza_file, qza_dir)

    organizing_files()

if __name__=="__main__":
    main()
