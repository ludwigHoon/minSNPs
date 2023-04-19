Tools used: 
- Dark-matter: https://antigenic-cartography.org/terry/consensus-making.html
- https://github.com/acorg/dark-matter/blob/master/bin/make-consensus.py

Deepsim needs data to be simulated to be in fasta format

```bash
conda activate bowtie
```

```python
from glob import glob
import subprocess

files = glob("./*.bam")

for file in files:
    id = file.split(".unmapped.bam")[0]
    id = id.replace("./", "")
    bashCommand = f"make-consensus.py --reference ../Mu50.fasta --bam {file} --id {id}"
    print(bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=open(f"{id}.fasta", "w+"))
    output, error = process.communicate()
    print(f"Done: {str(files.index(file)+1)}/{str(len(files))}")
```

From fasta to fast5

```python
from glob import glob
import subprocess

files = glob("./*.fasta")
#files = files[files.index("./ERR737075.fasta"):]

for file in files:
    #id = file.split(".fasta")[0]
    #id = id.replace("./", "")
    bashCommand = f"/home/lhoon/DeepSimulator/deep_simulator.sh -i {file}"
    print(bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=open(f"deep_sim.log", "a+"))
    output, error = process.communicate()
    print(f"Done: {str(files.index(file)+1)}/{str(len(files))}")

```

/usr/local/guppy_6.0.6/ont-guppy/bin/guppy_basecaller  -r --input_path ./fast5 \
    --save_path ./fastq2 -c dna_r9.4.1_450bps_hac.cfg \
    -x auto

bowtie2 -x mu50 -U fastq_runid_c2d19c211888bc09d8e077df271f325c911c1010_0_0.fastq --no-unal -p 1



readLines(pipe("sed -n '2~4p' STR1120-201113_DeepSimu/fastq/fail/fastq_runid_c2d19c211888bc
09d8e077df271f325c911c1010_0_0.fastq"), n = -1)

## Gradio, streamlit