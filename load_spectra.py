import tqdm
import os
import struct


def load_whole_pf2(path):
    '''
    Return a dict for experimental spectrum from .pf2 file
    '''

    f = open(path, "rb")
    spec_title = os.path.basename(path)[:-10]
    nSpec, lenTitle = struct.unpack("2i", f.read(8))
    pf2title = struct.unpack("%ds" %lenTitle, f.read(lenTitle))
    
    mpSpec = {}    
    for _ in tqdm.tqdm(range(nSpec)):
        scan_no, = struct.unpack("i",f.read(4))
        nPeak, = struct.unpack("i", f.read(4))
        peaks = []
        mz_int = struct.unpack(str(nPeak*2)+"d", f.read(nPeak*2*8))
        for i_peak in range(nPeak):
            mz = mz_int[i_peak*2]
            inten = mz_int[i_peak*2 + 1]
            peaks.append( (mz, inten) )
        
        #assert peaks == sorted(peaks, key=lambda x:x[0]) # pf2默认已经排过序了。
        max_inten=0
        if len(peaks)!=0:
            max_inten=max(peaks,key=lambda x:x[1])[1]
        
        nMix, = struct.unpack("i", f.read(4))
        nMaxCharge = 0
        for i_mix in range(nMix):
            precursor, = struct.unpack("d", f.read(8))
            nCharge, = struct.unpack("i", f.read(4))
            if nCharge > nMaxCharge: nMaxCharge = nCharge
            specname = f"{spec_title}.{scan_no}.{scan_no}.{nCharge}.{i_mix}.dta"
            spec_info = (specname, nCharge, precursor, max_inten)
            mpSpec[specname] = [spec_info, peaks]
    f.close()
    
    return pf2title, mpSpec

def pf2_loader(path):
    '''
    Return generator for experimental spectrum from .pf2 file
    '''
    
    with open(path, "rb") as f:
        spec_title = os.path.basename(path)[:-10]
        nSpec, lenTitle = struct.unpack("2i", f.read(8))
        pf2title = struct.unpack("%ds" %lenTitle, f.read(lenTitle))
        
        #mpSpec = {}    
        for _ in tqdm.tqdm(range(nSpec)):
            scan_no, = struct.unpack("i",f.read(4))
            nPeak, = struct.unpack("i", f.read(4))
            peaks = []
            mz_int = struct.unpack(str(nPeak*2)+"d", f.read(nPeak*2*8))
            for i_peak in range(nPeak):
                mz = mz_int[i_peak*2]
                inten = mz_int[i_peak*2 + 1]
                peaks.append( (mz, inten) )
            
            #assert peaks == sorted(peaks, key=lambda x:x[0]) #pf2默认已经排过序了。
            max_inten=0
            if len(peaks)!=0:
                max_inten=max(peaks,key=lambda x:x[1])[1]
            
            nMix, = struct.unpack("i", f.read(4))
            nMaxCharge = 0
            for i_mix in range(nMix):
                precursor, = struct.unpack("d", f.read(8))
                nCharge, = struct.unpack("i", f.read(4))
                if nCharge > nMaxCharge: nMaxCharge = nCharge
                specname = f"{spec_title}.{scan_no}.{scan_no}.{nCharge}.{i_mix}.dta"
                spec_info = (specname, nCharge, precursor, max_inten)
                yield specname, spec_info, peaks
    

def load_whole_mgf(path):
    '''
    Return a dict for experimental spectrum from .mgf file
    '''

    f = open(path, "r")
    
    mpSpec = {}    
    while "BEGIN IONS" in f.readline():    
        title = f.readline().split('=')[1].split('\n')[0]
        charge = f.readline().split('=')[1].split('\n')[0]
        rt = f.readline().split('=')[1].split('\n')[0]
        pepmass = f.readline().split('=')[1].split('\n')[0]

        peaks = ()
        line = f.readline()
        while not 'END IONS' in line:
            line = line.split('\n')[0].split(' ')
            peaks.append((float(line[0]), float(line[1])))
            line = f.readline()
        assert peaks == sorted(peaks, key=lambda x:x[0])
        spec_info=(title, charge, rt, pepmass)
        mpSpec[title] = [spec_info, peaks]

    f.close()
    
    return mpSpec
    