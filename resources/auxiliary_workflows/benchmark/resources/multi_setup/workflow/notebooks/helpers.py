method_palette = {"CliqueSNV": "#9B5DE5",
                  "cliquesnv": "#9B5DE5",
                  "HaploClique": "#FEE440", 
                  "haploclique": "#FEE440",
                  "HaploConduct": "#F15BB5", 
                  "haploconduct":"#F15BB5", 
                  "PredictHaplo":  "#00BBF9", 
                  "predicthaplo":  "#00BBF9", 
                  "ground_truth": "grey",
                  "ground truth": "grey",
                  }

method_order = ["CliqueSNV","HaploClique", "HaploConduct", "PredictHaplo"]

                  
def f_method_name(row):
    if row["method"] == 'cliquesnv':
        return "CliqueSNV"
    if row["method"] in ["haploclique"]:
        return "HaploClique"
    if row["method"] == 'haploconduct':
        return "HaploConduct"
    if row["method"] == 'predicthaplo':
        return "PredictHaplo"
    if row["method"] == 'ground_truth':
        return "ground truth"
    else: 
        print(row["method"])
