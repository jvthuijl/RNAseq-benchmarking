while read line; do 
    probename=$(echo $line | awk 'BEGIN {FS=";"}{print $1}' | grep -o -E '[0-9]+')
    rest=$(echo $line | awk 'BEGIN {FS=";"}{print $2}' | grep -o -E '[0-9]+')
    for i in $(echo $rest | awk 'BEGIN {FS="|"; OFS="\n"}{print $0}'); do
        echo $probename" ENST"$i
    done
done < $1
