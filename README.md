# Benchmarking RNA-seq tools
## addTranscriptsToProbe
The programs in this folder are used to add corresponding transcripts to unnamed probes from kallisto results. The file kallisto_probe_order.txt contains in which order the probe ids are sorted in kallisto. The file probe_to_transcript.out contains the probe id's and their corresponding transcripts. Both files are used in the program addTranscriptsToProbeOrder.sh, which returns a file in the same order as the kallisto results, with the probe names and corresponding transcripts.

## FixStrandInfo
This script fixes the strandinfo in the IlluminaProbe gtf file, as this
information was missing from the original gtf. The strand information is retrieved from a fasta file which contains the strand info for the probenames 

usage: parameter 1 is fileToRead = IlluminaMatcheRNAseqPRobes, paramter 2 is the fasta containing strand information.

## GenerateGraphs
Contains some of the R scripts used to create graphs for the results. Most importantly is the function CalcAbundance, which reads kallisto results in a folder, and combines them in a mean tpm and standard deviation for each transcript.

## SimulateReads
Contains the scripts used to simulate the groups base and diff with Polyester.

## Mutations
This program can be used to generate random snp's in a fasta file. It can be compiled on Linux and Unix systems using 

`g++ -std=c++14 -Wall main.cc -o mutations`

It can then be used as follows:

`./mutations {FASTAFileTOMutate}`

The output is caputered in {FASTAFileTOMutate}_mutated.fa

## removeFromFASTA
This program removes certain entities and their code from a FASTA file. It expects two argumtens. The first parameter is the FASTA
file itself, while the second parameter is the file containing the identifiers to remove.

 It can be compiled on Linux and Unix systems using:
 
`g++ -std=c++14 -Wall main.cc -o remFromFasta`

It can then be used as follows:

`./remFromFasta {FASTAFileTOMutate} {identifiers} > {outputFASTA}`

# Where is the rest?
Most of the other programs and data files are located in `/groups/umcg-wijmenga/tmp04/umcg-jvanthuijl`


