def write_mgf(wf, spec_info, peaks):    
    '''
    Write a spectrum with @specinfo as header and @peaks as array into @wf
    '''

    wf.write("BEGIN IONS\n")
    wf.write(f"TITLE={spec_info['TITLE']}\n")
    wf.write(f"CHARGE={spec_info['CHARGE']}\n")
    wf.write(f"RTINSECONDS={spec_info['RTINSECONDS']}\n")
    wf.write(f"PEPMASS={spec_info['PEPMASS']}\n")
    wf.write(f"SEQ={spec_info['SEQ']}\n")
    wf.write(peaks)
    wf.write("END IONS\n\n")
