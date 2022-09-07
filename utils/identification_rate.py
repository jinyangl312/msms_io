from ..load_spectra import get_mgf_titles
from ..load_results import get_seq_mod_from_pL_res
import pandas as pd

def calculate_scan_identification_rate(mgf_path, csv_path_list):
    """Calculate identification rate"""
    
    total_scans = set([x.split(".")[1] for x in get_mgf_titles(mgf_path)])

    total_valid_spec = pd.DataFrame()
    for csv in csv_path_list:
        spec = get_seq_mod_from_pL_res(csv)
        total_valid_spec = total_valid_spec.append(spec)

    total_valid_scans = set([x.split(".")[1] for x in total_valid_spec.index])

    return len(total_valid_scans)/len(total_scans)
