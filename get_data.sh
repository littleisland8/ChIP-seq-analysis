#get fasta 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O resources/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd resources

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

#get refFlat from UCSC
whet https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz -O resources/refFlat.txt

#get blacklist regions from ENCODE
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
