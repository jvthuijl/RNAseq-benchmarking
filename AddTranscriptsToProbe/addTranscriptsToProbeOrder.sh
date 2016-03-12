#!/bin/bash
probe_id="probe"
i=1;
for gene_id in $(cat kallisto_probe_order.txt); do
    MatchingTranscripts=($(awk "BEGIN {}{if (\$1==$gene_id){ print \$2}}" probe_to_transcript.out))
    if [ ${#MatchingTranscripts[@]} -eq 0 ]; then
        echo $probe_id
    else
        for transcript in ${MatchingTranscripts[@]}; do
            echo "$probe_id	$transcript"
        done 
    fi
    probe_id="probe_"$i 
    let i=$i+1
done;
