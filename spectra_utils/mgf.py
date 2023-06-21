import tqdm
import os
import re
from pathlib import Path
import numpy as np
from peaks_scripts.sequence_set import SequentialSet


def mgf_loader_unit(path, transform_peaks=True):
    '''
    Return a generator for spectra from .mgf text file @path
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    with open(path, "r") as f:
        while True:
            # Go to the next "BEGIN IONS"
            line = f.readline()
            if not line:  # EOF
                break
            while not "BEGIN IONS" in line:
                line = f.readline()
                if not line:  # EOF
                    break
            if not line:
                break

            # Parse spec headers
            line = f.readline()
            spec_info = dict()
            while "=" in line:
                line = re.split("=|\n", line)
                spec_info[line[0]] = line[1]
                line = f.readline()

            # Parse mz and intensity arrays
            if transform_peaks:
                peaks = []
                while not 'END IONS' in line:
                    line = re.split("\s|\n", line)
                    peaks.append((float(line[0]), float(line[1])))
                    line = f.readline()
            else:
                peaks = ""
                while not 'END IONS' in line:
                    peaks += line
                    line = f.readline()
            yield spec_info, peaks


def mgf_loader(path_list, transform_peaks=True):
    '''
    Return a generator for spectra from all .mgf text file in @path_list
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    if not isinstance(path_list, list):
        path_list = Path(path_list)
        assert path_list.is_dir()
        path_list = path_list.glob('*.mgf')

    for path in tqdm.tqdm(path_list):
        for spec_info, peaks in mgf_loader_unit(path):
            yield spec_info, peaks


def load_whole_mgf_unit(path, transform_peaks=True):
    '''
    Return a dict for all spectra from .mgf text file @path
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    if os.path.exists(f"{path}{'_t' if transform_peaks else '_k'}.npy"):
        return np.load(f"{path}{'_t' if transform_peaks else '_k'}.npy", allow_pickle=True).item()

    mpSpec = {}
    for spec_info, peaks in mgf_loader_unit(path, transform_peaks):
        mpSpec[spec_info["TITLE"]] = (spec_info, peaks)
    np.save(f"{path}{'_t' if transform_peaks else '_k'}.npy", mpSpec)
    return mpSpec


def load_whole_mgf(path_list, transform_peaks=True):
    '''
    Return a dict for all spectra from all .mgf text file in @path_list
    If @transform_peaks=True, mz and intensity arrays will be converted to float automatically;     otherwise they will be kept as a string.
    '''

    if not isinstance(path_list, list):
        path_list = Path(path_list)
        assert path_list.is_dir()
        path_list = path_list.glob('*.mgf')

    mpSpec = {}
    for path in tqdm.tqdm(path_list):
        mpSpec.update(load_whole_mgf_unit(path, transform_peaks))
    return mpSpec


class get_mgf_headers:
    '''
    Return a list of all headers contained in mgf.
    Used to check precursor evidence for results.
    '''

    def __init__(self, path_list, relative_error=20e-6, mixed_spectra=True):

        path_list = Path(path_list)
        assert path_list.is_dir()

        file_path_list = path_list.glob('*.mgf')

        self.mpSpec = {}
        for path in tqdm.tqdm(file_path_list, desc='Loading'):
            mpSpec_tmp = load_whole_mgf_unit(path, transform_peaks=False)
            mpSpec_tmp = {k: v[0] for k, v in mpSpec_tmp.items()}
            self.mpSpec.update(mpSpec_tmp)

        # Build up precursor info dict
        self.scan_headers_dict = dict()
        self.scan_charge_headers_dict = dict()
        self.PSM_headers_dict = dict()
        if mixed_spectra:
            for header in self.mpSpec.values():
                self.PSM_headers_dict[header["TITLE"]] = [
                    float(header["PEPMASS"])]

                scan_id = re.search(
                    ".*?\.\d+\.\d+(?=\.\d+\.\d+\.dta)", header["TITLE"]).group()
                if scan_id in self.scan_headers_dict:
                    self.scan_headers_dict[scan_id].append(
                        float(header["PEPMASS"]))
                else:
                    self.scan_headers_dict[scan_id] = [
                        float(header["PEPMASS"])]

                scan_charge_id = re.search(
                    ".*?\.\d+\.\d+\.\d+(?=\.\d+\.dta)", header["TITLE"]).group()
                if scan_charge_id in self.scan_charge_headers_dict:
                    self.scan_charge_headers_dict[scan_charge_id].append(
                        float(header["PEPMASS"]))
                else:
                    self.scan_charge_headers_dict[scan_charge_id] = [
                        float(header["PEPMASS"])]
        else:
            for header in self.mpSpec.values():
                self.PSM_headers_dict[header["TITLE"]] = [
                    float(header["PEPMASS"])]

                scan_id = re.search(
                    ".*?\.\d+\.\d+(?=\.\d+\.dta)", header["TITLE"]).group()
                if scan_id in self.scan_headers_dict:
                    self.scan_headers_dict[scan_id].append(
                        float(header["PEPMASS"]))
                else:
                    self.scan_headers_dict[scan_id] = [
                        float(header["PEPMASS"])]

                scan_charge_id = re.search(
                    ".*?\.\d+\.\d+\.\d+(?=\.dta)", header["TITLE"]).group()
                if scan_charge_id in self.scan_charge_headers_dict:
                    self.scan_charge_headers_dict[scan_charge_id].append(
                        float(header["PEPMASS"]))
                else:
                    self.scan_charge_headers_dict[scan_charge_id] = [
                        float(header["PEPMASS"])]
        print('Total PSMs:', len(self.PSM_headers_dict))
        print('Total charged scans:', len(self.scan_charge_headers_dict))
        print('Total scans:', len(self.scan_headers_dict))

        for k, v in self.PSM_headers_dict.items():
            self.PSM_headers_dict[k] = SequentialSet(v)
            self.PSM_headers_dict[k].set_relative_error(relative_error)
        for k, v in self.scan_headers_dict.items():
            self.scan_headers_dict[k] = SequentialSet(v)
            self.scan_headers_dict[k].set_relative_error(relative_error)
        for k, v in self.scan_charge_headers_dict.items():
            self.scan_charge_headers_dict[k] = SequentialSet(v)
            self.scan_charge_headers_dict[k].set_relative_error(relative_error)

        return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        return
