
# Database Curation for KoSCAPEbio

The Database Curation is a companion command-line tool for the KoSCAPEbio QIIME 2 plugin, designed to facilitate the bulk downloading of genomes, extraction of 16S rRNA and hypervariable regions, and preparation of data compatible with KoSCAPEbio. It allows users to expand or customize their databases for specific microbial community analyses.

## Features

- **Bulk Downloading**: Automates downloading genomes based on specified criteria.
- **Region Extraction**: Extracts 16S rRNA genes and creates hypervariable regions.
- **QIIME2 Compatibility**: Converts output files into QIIME 2 artifact format (.qza).
- **Customization**: Offers options for tailoring the database building process.
- **Standalone Flexibility**: Operates as a standalone tool or integrates with the KoSCAPEbio suite.

## Workflow Overview

Below is a visual overview of the workflow for database curation:
<div align="center">
    <img src="https://github.com/amar-cosic/q2-KoSCAPEbio/blob/master/q2_koscapebio/figures/db_workflow_overview.png" alt="Workflow Overview" width="600"/>
</div>

## Prerequisites

Before using this script, ensure Python 3.8 or higher is installed on your system. Additionally, this script depends on several third-party libraries not included in the Python Standard Library.

### Required Libraries:

All required libraries can either be installed via pip using the command: `pip install biopython pandas numpy scikit-bio` or by running: `pip install -r requirements.txt`. The `requirements.txt` file includes:

- `biom-format>=2.1.10`
- `pandas==1.5.3`
- `matplotlib>=3.5.0`
- `seaborn>=0.11.1`
- `scipy==1.10.0`
- `numpy>=1.21.0`
- `scikit-learn>=0.24.2`
- `Biopython==1.83`
- `scikit-bio==0.5.8`

## QIIME 2 Prerequisites

For utilizing the `-c --convert_qza` option, which converts output files into QIIME 2 artifact format (.qza), ensure [QIIME 2](https://qiime2.org) is installed and properly set up in your environment. Refer to the [official QIIME 2 documentation](https://docs.qiime2.org/) for installation instructions and environment setup.

## System Requirements

This tool is designed for Linux and macOS environments. Windows user may use WSL or VM boxes. 

## Quick Start

Before exploring all the customizable options available with `database_curation.py`, here's a quick start example to get you running with the basic functionality:

```bash
python3 database_curation.py --output ./custom_db_output --genome_number 5 --email your_email@example.com --convert_qza
```
This command initiates the download of 5 genomes, processes them to extract 16S rRNA and hypervariable regions, and outputs the data in QIIME 2 compatible format, all stored within the `./custom_db_output` directory.

## Installation

Clone the repository or download the `database_curation.py` script into your desired directory. In the same folder where the script is saved, store bla_database folder. This will ensure possibility to run bla-oxy checks.

## Detailed Parameter Descriptions

Below are the options you can use with `database_curation.py` to customize your database building process split between different categories.

#### Required
```
`-m`, `--list_genomes [LIST_GENOMES ...]`: Defines a list of specific genome identifiers or names to download and process, allowing for targeted database customization.
`-g`, `--genome_number`: Specifies the number of genomes to download randomly if `-l` is not specified. Useful for bulk processing.
`-o`, `--output`: Specifies the directory where output files will be saved. Default is the current directory.
`-e`, `--email EMAIL`: Your email address, required for accessing NCBI's Entrez databases via Biopython.
```
#### Recommended
```
`-b`, `--blaoxy_check`: Performs an additional verification step for Blaoxy.
`-l`, `--logging`: Enables logging for additional output and debugging information.
`-a`, `--api_key API_KEY`: API key for NCBI's Entrez databases when accessed via Biopython. Optional but recommended, especially if downloading a large number of genomes.
`-k`, `--consensus`: Generates a consensus sequence from the extracted 16S rRNA sequences, aiding in further analysis.
`--refseq`: Download only RefSeq assemblies if set.
```
#### Optional
```
`-d`, `--download_skip`: Skips the genome downloading step. Useful if you already have the genomes and need to modify other parameters.
`-c`, `--convert_qza`: Converts the final output into Qiime artifacts (.qza), making them directly usable in KoSCAPEbio.
`-s`, `--extraction_skip`: Skips the extraction of 16S and hypervariable regions. Useful when only file conversion to `.qza` format is needed.
`-2`, `--pairwise2`: Uses the deprecated `pairwise2` algorithm for alignments. Not recommended due to deprecation warnings but available for compatibility.
`-p`, `--primers_off`: Uses hardcoded nucleotide locations for the search of hypervariable regions instead of primers.
```
#### Standalone
```
`-t`, `--standalone`: Runs the tool in standalone mode without the KoSCAPEbio plugin.
`-x`, `--bla_oxy BLA_OXY`: Specifies the path for the `bla_oxy` database. Only use this in standalone mode if the `bla_database` is located outside the current directory.
```

Standalone parameters are only for usage outside of qiime2 and koscapebio plugin. 

## Example Usage

Below are example commands to help you get started with `database_curation.py` in both KoSCAPEbio-integrated and Standalone modes.

### Example 1: KoSCAPEbio-Integrated Mode (QIIME 2 Compatibility)
This command downloads 10 genomes for three specified Klebsiella species, extracts the 16S rRNA region, and converts the output to QIIME 2’s .qza format for easy integration with KoSCAPEbio:
```
python3 database_curation.py \
  --output ./database_output \
  --list_genomes "Klebsiella oxytoca" "Klebsiella grimontii" "Klebsiella variicola" \
  --genome_number 10 \
  --email your_email@example.com \
  --api_key YOURAPIKEY \
  --refseq \
  --convert_qza \
  --consensus
```
**Explanation**
- `--output ./database_output`:Specifies where to save output files.
- `--list_genomes`:Provides a list of specific species names to include in the database.
- `--genome_number 10`:Downloads 10 genomes (if available) for each specified species.
- `--refseq`:Downloads only genomes that are avalible in the refseq database
- `--convert_qza`:Converts the output to QIIME 2 `.qza` format for KoSCAPEbio.
- `--consensus`:Generates a consensus sequence for the 16S rRNA regions.

### Example 2: Standalone Mode
This command runs the tool independently (without KoSCAPEbio QIIME2 plugin), downloading the same genomes and saving the output in .fasta format. It specifies an external bla_oxy database for additional checks.
```
python3 database_curation.py \
  --output ./database_output \
  --list_genomes "Klebsiella oxytoca" "Klebsiella grimontii" "Klebsiella variicola" \
  --genome_number 10 \
  --email your_email@example.com \
  --api_key YOURAPIKEY \
  --refseq \  
  --convert_qza \
  --consensus \
  --standalone \
  --bla_oxy path_to_bla_db
```
**Explanation**
- `--standalone`:Runs the tool without the KoSCAPEbio QIIME 2 integration.
- `bla_oxy path_to_bla_db`:Specifies an external path to the `bla_oxy` database for bla_oxy checks.
- Other parameters (` --output`,` --list_genomes`,` --genome_number`, etc.) function as described in the integrated mode example.

## Output Explanation

The tool generates a structured output directory containing all processed data and intermediate files. When specifying an output directory, you’ll find three main folders within it:

### 1. Genomes
This folder stores the initial genomes downloaded as `.gbff` files. Inside the `Genomes` folder, there are two subfolders:

- **blast_results**: Contains `.out` files with the results of the `bla_oxy` check. If a genome fails the `bla_oxy` check, its `.out` file name will include `REMOVED`, indicating which genomes did not pass the check.
- **fasta**: Contains `.fasta` files, which are the `.gbff` genome files converted to FASTA format.

### 2. Output_files
The `Output_files` folder holds all main outputs generated by the tool, organized into specific subfolders:

- **FASTA**: This folder includes all output files in FASTA format, such as saved 16S sequences and removed sequences if consensus sequence building was selected. You’ll find files for each hypervariable region pair:
  - V3-V4, V4, and V3-V5 regions.
  - Each pair has two file versions, labeled as `"clean"` or `"duplicated"`, which reflect whether identical sequences were removed.
    - **Clean**: Use this version if the database is built from a single species.
    - **Duplicated**: Use this version when building the database from multiple species.

- **QZA**: Contains the FASTA files in QZA format, making them ready for direct use with QIIME 2.

- **Alignment file**: A `.txt` file from the in-silico PCR alignment process, containing detailed alignment information used in generating hypervariable regions.

### 3. PipelineResultsArchive
The `PipelineResultsArchive` folder serves as a logging and debugging repository. It includes two main subfolders:

- **CSV**: Contains summary tables in CSV format for tracking genome acquisition and 16S extraction results:
  - **genome_aquiring.csv**: This file has headers ``Specie Name`, `File Name`, `File Definition`, `BlaOxy Check`, `Percent Identity`, and `Saved?``, summarizing each genome’s details and `bla_oxy` results.
  - **16S_extraction.csv**: Summarizes the 16S extraction process, showing details such as ``File Name`, `File Definition`, `Number of 16S Extracted`, `Reason for No Extraction`, `Number of Keep Sequences`, `Number of Outlier Sequences`, and `Number of Different Sequences Saved``.

- **Logs**: Stores all generated log files, capturing steps like 16S extraction, file conversions, MAFFT alignments, and other key stages in the pipeline.

This folder structure is designed for easy navigation, enabling users to track and validate every step of the database curation process.

## Output Directory Structure

When you specify an output directory, the tool generates the following folder structure:

```
output_directory/
├── Genomes
│   ├── blast_results
│   │   ├── genome1_REMOVED.out
│   │   └── genome2.out
│   ├── fasta
│   │   ├── genome1.fasta
│   │   └── genome2.fasta
│   ├── genome1.gbff
│   └── genome2.gbff
├── Output_files
│   ├── FASTA
│   │   ├── 16S_sequences_clean.fasta
│   │   ├── 16S_sequences_duplicated.fasta
│   │   ├── V3-V4_clean.fasta
│   │   └── V3-V4_duplicated.fasta
│   ├── QZA
│   │   ├── 16S_sequences_clean.qza
│   │   ├── 16S_sequences_duplicated.qza
│   │   ├── V3-V4_clean.qza
│   │   └── V3-V4_duplicated.qza
│   └── alignment_results.txt
├── PipelineResultsArchive
│   ├── CSV
│   │   ├── genome_aquiring.csv
│   │   └── 16S_extraction.csv
│   └── Logs
│       ├── extraction.log
│       ├── conversion.log
│       └── alignment.log
```

## Troubleshooting

Here are some common issues that users may encounter, along with tips for resolving them:

- **API Rate Limits**:  
  When downloading genomes from NCBI, API rate limits may restrict the number of requests you can make in a short period. If you’re experiencing slow or interrupted downloads, try adding an API key using the `-a` or `--api_key` parameter. This can increase your request rate and make downloads more reliable. You can obtain an API key from NCBI by creating an account and generating a key [here](https://www.ncbi.nlm.nih.gov/account/).

- **Access Denied by NCBI**:  
  If you are not using an API key, NCBI may sometimes refuse access entirely, especially during high-traffic periods. Providing an API key is strongly recommended for consistent access and faster download speeds.

If you encounter other issues or need further support, please refer to the **Contact** section below.


## License

This tool is BSD licensed, as found in the LICENSE file.

## Citation

If you use Database curation script in your research, please cite it as follows:
`[Citation Here]`

## Contact

For help and support, please contact:

- **Name**: Amar Cosic
- **Email**: [amar.cosic995@gmail.com](mailto:amar.cosic995@gmail.com)
