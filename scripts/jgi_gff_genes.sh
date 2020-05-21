
perl -i -p -e 's/ID=gene_\d+;Name=jgi\.p\|([^\|]+)\|.+transcriptId=(\d+)/ID=$1|$1_$2/'  $1
