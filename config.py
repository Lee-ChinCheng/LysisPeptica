#from pathlib import Path

def get_config():
    return {
        "batch_size": 128,
        "num_epochs": 200,
        "lr": 4e-4,
        "seq_len": 52,
        "d_model": 160,
    }


def get_config_original_version():
    return {
        "batch_size": 128,
        "num_epochs": 200,
        "lr": 4e-4,
        "seq_len": 52,
        "d_model": 160,
        "datasource": 'dzjxzyd/UniRef50_len_0_50',
        "sequence_column": 'Reference sequence',
        "model_folder": "weights",
        "model_basename": "tmodel_",
        "preload": "latest",
        "tokenizer_file": "tokenizer.json",
    }

