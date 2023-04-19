# Pbsim command used
```
pbsim --hmm_model ~/pbsim2/data/R95.model --difference-ratio 23:31:46 --accuracy-mean 0.9 combined.fasta |&tee pbsim.log
```

# 4 line per reads in fastq, 1st line is read id, 2nd line is sequence, 3rd line is +, 4th line is quality score
```
for a in `ls *.fastq`
do
split -l 4000 $a ${a/.fastq/}"_" --additional-suffix=.fastq
done
```