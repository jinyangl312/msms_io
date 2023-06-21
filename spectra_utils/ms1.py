import tqdm
import os
import re
from pathlib import Path
import numpy as np
from peaks_scripts.sequence_set import SequentialSet


def ms1_loader_unit(path, transform_peaks=True):
    '''
    Return a generator from .ms1 text file @path
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    with open(path, "r") as f:
        line = f.readline()
        # Go to the next line starts with S
        if not line:  # EOF
            return
        while not line[0] == "S":
            line = f.readline()
            if not line:  # EOF
                break
        while True:
            if not line:
                return

            # Parse headers
            spec_info = dict()
            spec_info["scan_no"] = re.split("\t|\n", line)[1]
            line = f.readline()
            while line[0] == "I":
                line = re.split("\t|\n", line)
                spec_info[line[1]] = line[2]
                line = f.readline()

            # Parse mz and intensity arrays
            if transform_peaks:
                peaks = []
                for _ in range(int(spec_info["NumberOfPeaks"])):
                    line = re.split("\s|\n", line)
                    peaks.append((float(line[0]), float(line[1])))
                    line = f.readline()
            else:
                peaks = ""
                for _ in range(int(spec_info["NumberOfPeaks"])):
                    peaks += line
                    line = f.readline()
            yield spec_info, peaks


def load_whole_ms1_unit(path, transform_peaks=True):
    '''
    Return a dict for all .ms1 text file @path
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    if os.path.exists(f"{path}{'_t' if transform_peaks else '_k'}.npy"):
        return np.load(f"{path}{'_t' if transform_peaks else '_k'}.npy", allow_pickle=True).item()

    mpSpec = {}
    for spec_info, peaks in ms1_loader_unit(path, transform_peaks):
        mpSpec[spec_info["scan_no"]] = (spec_info, peaks)
    np.save(f"{path}{'_t' if transform_peaks else '_k'}.npy", mpSpec)
    return mpSpec


class load_whole_ms1:
    '''
    Return a list of all headers contained in mgf.
    TODO: Not sure if it was still used.
    '''

    def __init__(self, path_list, relative_error=20e-6, transform_peaks=True):

        if not isinstance(path_list, list):
            path_list = Path(path_list)
            assert path_list.is_dir()
            path_list = path_list.glob('*.ms1')

        self.mpSpec = {}
        for path in tqdm.tqdm(path_list, desc='Loading'):
            self.mpSpec[Path(path).stem] = load_whole_ms1_unit(
                path, transform_peaks)

        for _, ms1_unit_dict in tqdm.tqdm(self.mpSpec.items(), desc='Processing'):
            for scan_no, ms1_spec in ms1_unit_dict.items():
                ms1_unit_dict[scan_no] = SequentialSet(
                    np.array(ms1_spec[1])[:, 0])
                ms1_unit_dict[scan_no].set_relative_error(relative_error)
        return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        return
