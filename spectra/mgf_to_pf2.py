import tqdm
from msms_io.mass.AAMass import aamass
import re
import struct
from msms_io.spectra.mgf import mgf_loader_unit
import pathlib
import os


def transform_mgf_to_pf2(mgf_dir):
    '''
    Transfer all mgf files in @mgf_dir into corresponding pf2 and pf2idx files
    Used for pBuild visualiaztion
    '''
    if not isinstance(mgf_dir, list):
        mgf_dir = pathlib.Path(mgf_dir)
        assert mgf_dir.is_dir()

    for mgf_path in tqdm.tqdm(mgf_dir.glob('*.mgf')):
        transform_mgf_to_pf2_unit(str(mgf_path))


def transform_mgf_to_pf2_unit(mgf_path):
    '''
    Transfer @mgf_path into corresponding pf2 and pf2idx files
    Sub-function for pBuild visualiaztion
    '''
    assert '_' in pathlib.Path(mgf_path).stem
    pf2idx_path = mgf_path.replace('.mgf', '.pf2idx')
    pf2_path = mgf_path.replace('.mgf', '.pf2')
    pf2title = pathlib.Path('_'.join(mgf_path.split('_')[:-1])).stem

    if os.path.exists(pf2idx_path) and os.path.exists(pf2_path):
        print('skipping pf2')
        return

    mgf_loader = mgf_loader_unit(mgf_path, transform_peaks=True)

    # Group mgf by scan_id, merge mixed spectra
    mgf_scan_dict = dict()
    mgf_scan_no_list = list()
    for spec_info, peaks in tqdm.tqdm(mgf_loader):
        scan_no = int(re.search("\d+(?=\.\d+\.\d+\.\d+\.dta)",
                      spec_info['TITLE']).group())
        if scan_no not in mgf_scan_dict:
            mgf_scan_dict[scan_no] = [(
                scan_no, len(peaks),
                [y for x in peaks for y in x]
            )]
            mgf_scan_no_list.append(scan_no)
        mgf_scan_dict[scan_no].append((
            float(spec_info['PEPMASS']),
            int(re.search("\d+", spec_info['CHARGE']).group())
        ))

    # Write pf2 and pf2idx
    pf2idx_file = open(pf2idx_path, 'wb')
    pf2_file = open(pf2_path, 'wb')

    offset = 0
    pf2_file.write(struct.pack("2i", len(mgf_scan_no_list), len(pf2title)))
    pf2_file.write(struct.pack("%ds" %
                   len(pf2title), bytes(pf2title, 'utf-8')))
    offset += 2*4 + len(pf2title)*1

    for scan_no in tqdm.tqdm(mgf_scan_no_list):
        pf2idx_file.write(struct.pack("i", scan_no))  # scan_no
        pf2idx_file.write(struct.pack("i", offset))  # offset

        scan_info = mgf_scan_dict[scan_no]
        pf2_file.write(struct.pack("i", scan_info[0][0]))  # scan_no
        pf2_file.write(struct.pack("i", scan_info[0][1]))  # nPeak
        pf2_file.write(struct.pack(
            str(scan_info[0][1]*2)+"d", *scan_info[0][2]))  # peaks
        pf2_file.write(struct.pack("i", len(scan_info)-1))  # nMix
        offset += 3*4 + scan_info[0][1]*2*8
        for i in range(1, len(scan_info)):
            pf2_file.write(struct.pack("d", scan_info[i][0]))  # precursor
            pf2_file.write(struct.pack("i", scan_info[i][1]))  # nCharge
        offset += (len(scan_info)-1) * (4+8)

    pf2idx_file.close()
    pf2_file.close()
