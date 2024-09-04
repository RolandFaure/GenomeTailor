# GenomeTailor

GenomeTailor edits the assembly graph in GFA format to make sure all reads align end-to-end on the graph. It also deletes the regions that are not covered by any reads. Currently under development and not very stable, if you need to perfect your assemblies I advise you write to me.

## Installation

### Dependencies

To compile the project you will need CMake>=3.8.12 and GCC (tested with GCC=11.3.1).
You will need to have `racon`, `minimap2`, `raven` and optionnally `minigraph`.

### Download & install

```
git clone https://github.com/RolandFaure/GenomeTailor.git
cd GenomeTailor
mkdir build && cd build
cmake ..
make
```

Test the installation on the test instance:
```
cd test
../build/GenomeTailor -i assembly.gfa -r mock_reads.fasta -o test.gfa
```

## Usage

```
build/GenomeTailor 
SYNOPSIS
        build/GenomeTailor -i <input_assembly> -r <input_reads> -m <mode> [-e <output_errors>] [-o
                           <output_assembly>] [-d <output_non_duplexed_reads>] [-p
                           <path-to-tmp-folder>] [-b <minimum-number-of-reads>] [-t <threads>] [-g
                           <gaf_file>] [--minigraph <minigraph>] [--minimap2 <minimap2>] [--racon
                           <racon>] [--path-to-raven <path-to-raven>] [-h] [-v]

OPTIONS
        -i, --input_assembly
                    input assembly in gfa format

        -r, --input_reads
                    input reads in fasta/q format

        -m, --mode  mode: correct or detect
        -e, --output_errors
                    output file describing the errors found in the assembly

        -o, --output_assembly
                    output assembly in gfa format (required if correct mode)

        -d, --output_non_duplexed_reads
                    file to output a file of non-duplexed reads

        -p, --path-to-tmp-folder
                    path to a temporary folder where the intermediate files will be stored [./]

        -b, --minimum-number-of-reads
                    minimum number of reads to support a breakpoint [5]

        -t, --threads
                    number of threads to use for minigraph [1]

        -g, --gaf_file
                    gaf file if already computed (NO SECONDARY ALIGNMENTS). Will be generated with
                    minigraph if not provided

        --minigraph path to minigraph
        --minimap2  path to minimap2
        --racon     path to racon
        --path-to-raven
                    path to raven

        -h, --help  print this help message and exit
        -v, --version
                    print version information and exit

```
