import csv

#>seq0|Conc=54.61
#ILPWKCPWWPCRR

input='Hemo_predi/H10_8DL49/test1_173.fa'
output='LysisPeptica/testset/test1.fa'
raw_data='LysisPeptica/raw_data/ulin2985.csv'

test_id=[]
with open(input, 'r') as f:
    for l in f:
        l=l.strip()
        #>16939	12.0_116.899_62.54	0.487
        #FLPIAGKLLSGLSGLL
        if l[0]=='>':
            rsid=l.split('\t')[0][1:]
            if rsid[0]=='-':
                sid=rsid[1:]
            else:
                sid=rsid        
            test_id.append(sid)
test_id=set(test_id)

with open(output, 'w') as op:

    with open(raw_data, newline='') as csvfile:
        reader = csv.reader(csvfile)
        #id,sequence,lysis_group,lysis_value,concentration,concentration_uM
        for row in reader:
            sid = int(row[0])       
            if sid in test_id:
                #print(row)
                content=f'>{sid}|Conc={row[4]}\tConc(uM)={row[5]}'

