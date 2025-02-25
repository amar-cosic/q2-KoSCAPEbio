import cmd
from q2_types.feature_data import FeatureData, Sequence
from qiime2.plugin import (Plugin, Str, Properties, Int, Float, Range,
                           Bool, Choices, Metadata)
from qiime2.plugin import SemanticType
from qiime2.plugin.model import File, DirectoryFormat
import sys

if sys.version_info >= (3, 9):
    import importlib.resources as resources
else:
    import importlib_resources as resources

import q2_koscapebio
from q2_types.feature_data import FeatureData, Sequence, DNAFASTAFormat, DNAIterator
import subprocess
from q2_types.feature_table import FeatureTable, Frequency
from qiime2 import Metadata
import qiime2
import biom
from biom import Table
from q2_types.feature_data import DNASequencesDirectoryFormat
import os
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import tempfile
import gzip
import shutil
import logging
import seaborn as sns


plugin = qiime2.plugin.Plugin(
    name='koscapebio',
    version=q2_koscapebio.__version__,
    website="https://github.com/amar-cosic/q2-KoSCAPEbio",
    package='q2_koscapebio',
    short_description='KoSCAPEbio: Klebsiella  oxytoca Species Complex Analysis and Profiling Enabler',
    description=('KoSCAPEbio is a QIIME 2 plugin designed for the exploration and '
                 'analysis of the Klebsiella oxytoca species complex (KoSC) presence within '
                 'microbial communities.')
)

def presence_check(rep_seq_path: str, 
                table_path: str, 
                perc_identity: float = 1,
                temp_dir: Str = None,
                user_db: Str = None,
                strand: Str = "plus") -> (biom.Table, DNASequencesDirectoryFormat, DNASequencesDirectoryFormat):

    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    
    if temp_dir is not None:
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
    if user_db is not None:
        _, ext = os.path.splitext(user_db)
        if ext.lower() == ".fasta" or ext.lower() == ".gz":
            if ext.lower() == ".gz":
                base_name, base_ext = os.path.splitext(os.path.splitext(user_db)[0])
                if base_ext.lower() != ".fasta":
                    raise ValueError("Invalid file format for user_db. Please provide a .fasta, .fasta.gz, or .qza file.")
                unzipped_path = os.path.join(temp_dir, "user_db.fasta")
                with gzip.open(user_db, "rb") as f_in:
                    with open(unzipped_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                db_path = os.path.join(temp_dir, "user_db.qza")
                cmd = [
                    "qiime", "tools", "import",
                    "--input-path", unzipped_path,
                    "--output-path", db_path,
                    "--type", "FeatureData[Sequence]"
                ]
                subprocess.run(cmd, check=True)
                os.remove(unzipped_path)
            else:
                db_path = os.path.join(temp_dir, "user_db.qza")
                cmd = [
                    "qiime", "tools", "import",
                    "--input-path", user_db,
                    "--output-path", db_path,
                    "--type", "FeatureData[Sequence]"
                ]
                subprocess.run(cmd, check=True)

        elif ext.lower() == ".qza":
            db_path = user_db

        else:
            raise ValueError("Invalid file format for user_db. Please provide a .fasta, .fasta.gz, or .qza file.")

    else:
        db_path = resources.files("q2_koscapebio.data").joinpath("database.qza")

    clustered_table_path = os.path.join(temp_dir, "positive_table.qza")
    clustered_seq_path = os.path.join(temp_dir, "positive_seq.qza")
    unclustered_seq_path = os.path.join(temp_dir, "negative_seq.qza")

    cmd = [
        "qiime", "vsearch", "cluster-features-closed-reference",
        "--i-sequences", rep_seq_path,
        "--i-table", table_path,
        "--i-reference-sequences", db_path,
        "--p-perc-identity", str(perc_identity),
        "--p-strand", strand,
        "--o-clustered-table", clustered_table_path,
        "--o-clustered-sequences", clustered_seq_path,
        "--o-unmatched-sequences", unclustered_seq_path
    ]
    subprocess.run(cmd,check=True)

    clustered_table = qiime2.Artifact.load(clustered_table_path).view(biom.Table)
    clustered_seq = qiime2.Artifact.load(clustered_seq_path).view(DNASequencesDirectoryFormat)
    unclustered_seq = qiime2.Artifact.load(unclustered_seq_path).view(DNASequencesDirectoryFormat)

    return (clustered_table, clustered_seq, unclustered_seq)

plugin.methods.register_function(
    function=presence_check,
    inputs={},
    parameters={
        "rep_seq_path": Str,
        "table_path": Str,
        "perc_identity": Float,
        "temp_dir": Str,
        "user_db": Str,
        "strand": Str % qiime2.plugin.Choices(["plus","both"])
    },
    outputs=[
        ("clustered_table", FeatureTable[Frequency]),
        ("clustered_seq", FeatureData[Sequence]),
        ("unclustered_seq", FeatureData[Sequence])
    ],
    input_descriptions={},
    parameter_descriptions={
        "rep_seq_path": "Path for representative sequences in '.qza' format.",
        "table_path": "Path for the feature table in '.qza' format.",
        "perc_identity": "The percentage of identity for clustering. Between 0 and 1",
        "temp_dir":"(Optional) Directory for temporary files. The files will be deleted. If not provided, a temporary directory will be created.",
        "user_db":"User defined database",
        "strand":"Search forward strand or both the forward and reverse sequence. "
    },
    output_descriptions={
        "clustered_table": "Feature table with positive OTU/ASV in samples.",
        "clustered_seq": "Positive sequences that could be found in the samples.",
        "unclustered_seq": "The resulting unmatched sequences."
    },
    name="presence_check",
    description="PresenceCheck scans the input representative sequences and feature table against a custom database to identify the presence of Klebsiella oxytoca species complex (KoSC) members. Utilizing VSEARCH algorithms, it offers  detection based on specified percentage identity criteria."
)

def abundance_profile(raw_table: str, 
                qiime_table: str, 
                work_dir: str,) -> (biom.Table):
    
    raw_table_pth=f"{work_dir}/koscapebio_raw_table/"
    qiime_table_pth=f"{work_dir}/qiime_table/"
    cmd_raw= [
        "qiime","tools","export",
        "--input-path",raw_table,
        "--output-path",f"{raw_table_pth}"
    ]
    subprocess.run(cmd_raw,check=True)
    cmd_q= [
        "qiime","tools","export",
        "--input-path",qiime_table,
        "--output-path",f"{qiime_table_pth}"
    ]
    subprocess.run(cmd_q,check=True)

    cmd_r_c=[
        "biom",
        "convert",
        "-i",f"{raw_table_pth}feature-table.biom",
        "-o",f"{raw_table_pth}koscapebio_raw_table.tsv",
        "--to-tsv"
    ]
    subprocess.run(cmd_r_c,check=True)
    cmd_q_c=[
        "biom",
        "convert",
        "-i",f"{qiime_table_pth}/feature-table.biom",
        "-o",f"{qiime_table_pth}/qiime_table.tsv",
        "--to-tsv"
    ]
    subprocess.run(cmd_q_c,check=True)

    output_path=os.path.join(work_dir, "final_table.tsv")
    raw_table_tsv=f"{raw_table_pth}koscapebio_raw_table.tsv"
    qiime_table_tsv=f"{qiime_table_pth}qiime_table.tsv"
    try:
        total_abundance_df = pd.read_csv(qiime_table_tsv, sep='\t', header=1)
        total_abundance_df = total_abundance_df.drop(total_abundance_df.columns[0], axis=1)
        total_abundances = total_abundance_df.sum()

        raw_abundance_df = pd.read_csv(raw_table_tsv, sep='\t', header=1)
        feature_ids = raw_abundance_df.iloc[:, 0]  

        species_names = feature_ids.apply(lambda x: x if 'shared' in x.lower() else x.split('_')[1])

        raw_abundance_data = raw_abundance_df.drop(raw_abundance_df.columns[0], axis=1)
        raw_abundance_data['Species'] = species_names  
        grouped_abundance = raw_abundance_data.groupby('Species').sum()

        relative_abundance_data = grouped_abundance.div(total_abundances, axis=1)

        relative_abundance_data.to_csv(output_path, sep='\t')
        print(f"Relative abundance data saved to {output_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

    cmd_r=[
        "biom",
        "convert",
        "-i",f"{output_path}",
        "-o",f"{work_dir}/rel_table.biom",
        "--to-hdf5"
    ]
    subprocess.run(cmd_r,check=True)

    data = pd.read_csv(output_path, sep='\t', index_col=0)
    print(data)

    data = data.T

    plt.figure(figsize=(10, 8))
    sns.heatmap(data, yticklabels=True, xticklabels=data.columns,cmap='viridis',annot=True)
    plt.xlabel('Sample Names')
    plt.ylabel('Features')
    plt.title('Heatmap of Transposed TSV Data')
    plt.tight_layout()

    plt.savefig(f'{work_dir}/heatmap.png')

    temp_dir = tempfile.mkdtemp()
    temp_test_file_path = os.path.join(temp_dir, "test.qza")
    cmd= [
        "qiime","tools","import",
        "--input-path",f"{work_dir}/rel_table.biom",
        "--output-path",temp_test_file_path,
        "--type","FeatureTable[Frequency]"
    ]
    subprocess.run(cmd,check=True)
    temp_test_file = qiime2.Artifact.load(temp_test_file_path).view(biom.Table)
    return (temp_test_file)

plugin.methods.register_function(
    function=abundance_profile,
    inputs={},
    parameters={
        "raw_table": Str,
        "qiime_table": Str,
        "work_dir": Str,
    },
    outputs=[
        ("relative_abundance", FeatureTable[Frequency]),
    ],
    input_descriptions={},
    parameter_descriptions={
        "raw_table": "Path to the raw table output from the KoSCAPEbio analysis. This table should contain the initial analysis results generated by KoSCAPEbio.",
        "qiime_table": "Path to the feature table generated by QIIME 2. This table is the result of QIIME 2's processing pipeline and is used for comparative analysis with the KoSCAPEbio output.",
        "work_dir": "Directory path where intermediate files will be stored. These files are generated during the visualization process and may be useful for detailed analysis or debugging.",
    },
    output_descriptions={
        "relative_abundance": "The resulting feature table indicating the relative abundance of features. This table is saved as a QIIME 2 artifact (.qza) and is ready for further analysis or visualization within the QIIME 2 framework.",
    },
    name="abundance_profile",
    description=("AbundanceProfile calculates and visualizes the relative abundance of microbial species, integrating KoSCAPEbio findings with broader QIIME 2 microbial community analyses.")
)