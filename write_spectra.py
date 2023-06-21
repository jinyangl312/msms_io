from pathlib import Path
import tqdm
from theoretical_peaks.AAMass import aamass
import re
from .load_spectra import mgf_loader_unit
from .spectra_utils.mgf_to_pf2 import *


def write_mgf(writer, header, peaks):
    '''
    Write a spectrum with @header and @peaks into @writer
    Sub-function for mgf manipulation
    '''

    writer.write("BEGIN IONS\n")
    for key, item in header.items():
        writer.write(f"{key}={item}\n")
    if isinstance(peaks, str):
        writer.write(peaks)
    else:
        processed_scan = ["{:.5f} {:.1f}".format(
            line[0], line[1]) for line in peaks]
        writer.write('\n'.join(processed_scan))
        writer.write("\n")

    writer.write("END IONS\n")


def shift_precursor_masses_012(input_list_dir, output_dir):
    '''
    variate spectrum in mgf files in @input_list_dir by adding precursor mass up to 1 and 2 Da lighter
    Used to preprocess pXtract results before search with pLink
    '''

    input_list_dir = Path(input_list_dir)
    assert input_list_dir.is_dir()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for path in tqdm.tqdm(input_list_dir.glob('*.mgf')):
        with open(output_dir.joinpath(path.name), "w") as fout:
            for spec_info, peaks in mgf_loader_unit(path):
                spec_info["PEPMASS"] = float(spec_info["PEPMASS"])
                charge = int(spec_info["CHARGE"].strip('+'))
                for i in range(0, 3):
                    spec_info_new = spec_info.copy()
                    spec_info_new["TITLE"] = f'{spec_info_new["TITLE"].split(".dta")[0]}.{i}.dta'
                    spec_info_new["PEPMASS"] = "{}".format(
                        spec_info["PEPMASS"] - i * aamass.mass_isotope / charge)
                    write_mgf(fout, spec_info_new, peaks)


def add_scans_info(input_list_dir, output_dir):
    '''
    Add scan header info for xi search
    Used to preprocess pXtract results before search with xi
    '''

    input_list_dir = Path(input_list_dir)
    assert input_list_dir.is_dir()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for path in tqdm.tqdm(input_list_dir.glob('*.mgf')):
        with open(output_dir.joinpath(path.name), "w") as fout:
            for spec_info, peaks in mgf_loader_unit(path):
                spec_info["SCANS"] = re.search(
                    "\d+(?=\.\d+\.\d+\.dta)", spec_info["TITLE"]).group()
                write_mgf(fout, spec_info, peaks)
