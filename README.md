# LysisPeptica
A docker for therapeutic peptides, hemolysis risk prediction. The pipeline includes selectable hemolysis% threshold for prediction 5, 10, 20 or 30. Default is 10.

### Main Goal
Therapeutic Peptides or Antimicrobial peptides (AMPs) with high hemolysis risk are unsuitable as drugs, but experimental validation is costly. Our model enables early in silico screening to filter out hemolytic candidates.

Our web-based AI predictor, LysisPeptica is accessible at https://axp.iis.sinica.edu.tw/hemolysis/

---

### Step by step to deploy this docker

1. download this repo folder in your local machine

2. go to folder directory
~$ cd LysisPeptica

3. build docker image
~$ docker build -t lysispeptica .

(wait for some minutes)

4. test input and output
~$ docker run --rm \
  -v /your_path/your.fa:/app/data/input.fasta \
  -v /your_output_folder:/app/output \
  lysispeptica --thr_id {hemolysis%}

replace "/your_path/your.fa" by your input fasta path
replace "/your_output_folder" by your preferred folder for output csv file
{hemolysis%} wil be 5,10,20,30. default is 10