import os, argparse
import numpy as np, pandas as pd
#os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" #warning of CPU , GPU utilization
import tensorflow as tf
import torch
from tokenizers import Tokenizer
from .config import get_config
from .model import build_transformer
from .func import pc6_8d_encode, add_conc_on_pepbert_array
from .func import CustomModel, GlobalMinPooling1D



###==== pepBERT encoding ====
 
# 2) recall tokenizer.json and load the tokenizer
tokenizer_path = 'hf_hub/tokenizer.json'
tokenizer = Tokenizer.from_file(tokenizer_path)

# 3) recall model weights (has downloaded from hf)
weights_path = 'hf_hub/tmodel_16.pt'

# 4) Initialize the model structure and load the weights
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_built() or torch.backends.mps.is_available() else "cpu"
config = get_config()

model = build_transformer(
    src_vocab_size=tokenizer.get_vocab_size(),
    src_seq_len=config["seq_len"],
    d_model=config["d_model"]
)
state = torch.load(weights_path, weights_only=True, map_location=torch.device(device)) 
#default is weights_only=False, will raise a FutureWarning
model.load_state_dict(state["model_state_dict"])
model.eval()



def pbert_encode( seqli, tar_len):
    opli=[]
    for seq in seqli:       
        en_ids = ([tokenizer.token_to_id("[SOS]")] + tokenizer.encode(seq).ids + [tokenizer.token_to_id("[EOS]")])
        en_ids = en_ids + [tokenizer.token_to_id("[PAD]")] * (tar_len-len(seq))
        input_ids = torch.tensor([en_ids], dtype=torch.int64)
        with torch.no_grad():
            # Attention mask (1 for tokens, 0 for PADs)
            encoder_mask = (input_ids != tokenizer.token_to_id("[PAD]")).unsqueeze(1).unsqueeze(2).long()
            # Forward pass through the encoder to get token embeddings
            emb = model.encode(input_ids, encoder_mask)
            # Remove first ([SOS]) and last ([EOS]) embeddings
            emb_no_first_last = emb[:, 1:-1, :] 

            # Apply average pooling over the remaining tokens to get a fixed-size vector
            # emb_avg = emb_no_first_last.mean(dim=1) #->torch.Size([1, 320])

            # save in numpy
            ar_np = emb_no_first_last.squeeze(0).numpy() 
            opli.append(ar_np)

    # opli can't be numpy array, cause will add conc into 161 depth      
    return opli




def read_fasta(filepath: str, df_conc: int = 50):

    ids, concs, seqs = [], [], []
    seq_id = None
    seq_conc = None
    # ---------- Parse FASTA ----------
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
              
                # Parse new header
                header = line[1:]
                # Ex: seq1|Conc=14.3
                parts = header.split("|")
                seq_id = parts[0]
                ids.append(seq_id)
                p2 = parts[1].split("=")[1]
                try:
                    seq_conc = float(p2)
                except:      
                    seq_conc = 50 #assign default conc 50 ug/ml
                concs.append(seq_conc)
            else:
                seqs.append(line)

    return ids, seqs, concs


def read_fasta_slice(filepath: str, window: int = 49, df_conc: int = 50):
    """
    Read FASTA and output 3 lists:
      - ids: fragment IDs (e.g. seq1_1_10)
      - concs: float concentrations
      - seqs: fragment sequences
    """
    ids, concs, seqs = [], [], []
    seq_id = None
    seq_conc = None
    seq_lines = []
    # ---------- Parse FASTA ----------
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Yield previous record
                if seq_id is not None:
                    full_seq = "".join(seq_lines)
                    _slice_sequence(seq_id, seq_conc, full_seq, window, ids, concs, seqs)

                # Parse new header
                header = line[1:]
                # Example: seq1|Conc=14.3
                parts = header.split("|")
                seq_id = parts[0]
                p2 = parts[1].split("=")[1]

                try:
                    seq_conc = float(p2)
                except:      
                    seq_conc = 50 #assign default conc 50 ug/ml
                seq_lines = []
            else:
                seq_lines.append(line)

        # Final record
        if seq_id is not None:
            full_seq = "".join(seq_lines)
            _slice_sequence(seq_id, seq_conc, full_seq, window, ids, concs, seqs)

    return ids, seqs, concs


def _slice_sequence(seq_id, conc, full_seq, window, ids, concs, seqs):
    """
    Slice sequence using non-overlapping sliding windows.
    Append results to ids, concs, seqs lists.
    """
    seq_len = len(full_seq)

    # Case 1: sequence <= window → keep original
    if seq_len <= window:
        ids.append(seq_id)
        concs.append(conc)
        seqs.append(full_seq)
        return

    # Case 2: slice into fragments
    start = 1  # human-readable index
    while start <= seq_len:
        end = min(start + window - 1, seq_len)
        frag = full_seq[start-1:end]  # Python index is 0-based
        frag_id = f"{seq_id}_{start}_{end}"
        ids.append(frag_id)
        concs.append(conc)
        seqs.append(frag)
        start = end + 1



def encoded_policy(encoding, mdpath, seqli, ugmlli):
    if encoding =='pc6zs':
        X_input  = pc6_8d_encode( seqli, ugmlli, encoding, 49 ) 
        
    elif encoding =='pc6norm':
        X_input  = pc6_8d_encode( seqli, ugmlli, encoding, 49 ) 

    elif encoding =='pepbert_ugml':
        d160li = pbert_encode(seqli, 49)
        X_input = add_conc_on_pepbert_array(encoding, d160li, seqli, ugmlli)

    else: # encoding =='pepbert_um':
        d160li = pbert_encode(seqli, 49)
        X_input = add_conc_on_pepbert_array(encoding, d160li, seqli, ugmlli)
        

    model = tf.keras.models.load_model(mdpath,
                                        custom_objects={'CustomModel': CustomModel, 'tf': tf,
                                                        'GlobalMinPooling1D':GlobalMinPooling1D},compile=False)
    pred_probs = model.predict(X_input)
    #print(type(pred_probs)) #'numpy.ndarray', 2depth
    #[[0.509743   0.49025705]
    #[0.7774866  0.22251332]]
    pred_probs=pred_probs[:, 1] #ex: #[0.49025705 0.22251332]]
    return pred_probs 


def ensemble_prob(prob_li):
    prob_li = np.array(prob_li)
    ens_prob = prob_li.mean(axis=0)
    return(ens_prob)



def predict(thr_id: int, input_fasta: str, output_csv: str):
    #Runs prediction pipeline using the selected models, based on selected hemolysis threshold.
    # models and mapping encodings policy as below
    md_policy={
        5:[ ('md845613_5240chatt.keras','pepbert_um') ],

        10:[('p791_836_cnn_zs_5544.keras','pc6zs'), 
            ('p798_796_cnn2_zs_5545.keras','pc6zs'),
            ('p763_843_5950_3p1bn_ugml2std.keras','pepbert_ugml'),
            ('p843_750_5041chatt_ugml2std.keras','pepbert_ugml') ],

        20:[ ('md742811_1bnMLP.keras','pepbert_um'),('md749791_2p2bnMLP.keras','pepbert_um') ],
        30:[ ('md871776_3bnMLP.keras','pepbert_um') ]}

    #read input fasta
    #idli, seqli, ugmlli = read_fasta_slice(input_fasta, 25, 50)
    idli, seqli, ugmlli = read_fasta(input_fasta, 50)
   

    ens_li=[]
    if thr_id==5:
        for p in md_policy[thr_id]:
            mdpath = f'/app/models/thr5/{p[0]}'            
            pred_probs = encoded_policy(p[1], mdpath, seqli, ugmlli)
            ens_li.append(pred_probs)
        
    elif thr_id==10:
        for p in md_policy[thr_id]:
            mdpath = f'/app/models/thr10/{p[0]}'        
            pred_probs = encoded_policy(p[1], mdpath, seqli, ugmlli)
            ens_li.append(pred_probs)
        
    elif thr_id==20:
        for p in md_policy[thr_id]:
            mdpath = f'/app/models/thr20/{p[0]}'     
            pred_probs = encoded_policy(p[1], mdpath, seqli, ugmlli)
            ens_li.append(pred_probs)
        
    elif thr_id==30:
        for p in md_policy[thr_id]:
            mdpath = f'/app/models/thr30/{p[0]}'                      
            pred_probs = encoded_policy(p[1], mdpath, seqli, ugmlli)  
            ens_li.append(pred_probs) 
    
    else: #thr_id==10:
        for p in md_policy[thr_id]:
            mdpath = f'/app/models/thr10/{p[0]}'           
            pred_probs = encoded_policy(p[1], mdpath, seqli, ugmlli)
            ens_li.append(pred_probs)

    ens_li = np.array(ens_li)
    ens_probs = ens_li.mean(axis=0)  
    #print(ens_probs)


    thrli=[]
    for i in range(len(idli)):
        thrli.append(thr_id)
    binaryli=[]
    for i in ens_probs:
        if i >=0.5: binaryli.append('Yes')
        else:       binaryli.append('No')
   

    #print(f"Saving results to {output_csv}...")
    results_df = pd.DataFrame({
        'PEPTIDE': idli,
        'HEMOLYSIS THRESHOLD(%)': thrli,
        'DOSAGE(DEFAULT=50µg/ml)': ugmlli,
        'SCORE': ens_probs,
        'HEMOLYSIS': binaryli})
    
    results_df.to_csv(output_csv, index=False, float_format='%.4f')
  
       
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Peptide Hemolysis Prediction")
    parser.add_argument('--input', type=str, required=True, help="Path to the input FASTA file.")
    parser.add_argument('--output', type=str, required=True, help="Path for the output CSV file.")
    parser.add_argument(
        '--thr_id', 
        type=int, 
        default=10, 
        choices=(5, 10, 20, 30), 
        help="select hemolysis% threshold for prediction 5, 10, 20 or 30. Default is 10"
    )
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    predict(thr_id=args.thr_id, input_fasta=args.input, output_csv=args.output)



