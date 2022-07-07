def write_mgf(writer, header, peaks):    
    '''
    Write a spectrum with @header and @peaks into @wf
    '''

    writer.write("BEGIN IONS\n")
    writer.write(f"TITLE={header['TITLE']}\n")
    writer.write(f"CHARGE={header['CHARGE']}\n")
    writer.write(f"RTINSECONDS={header['RTINSECONDS']}\n")
    writer.write(f"PEPMASS={header['PEPMASS']}\n")
    writer.write(f"SEQ={header['SEQ']}\n")
    if isinstance(peaks, str):
        writer.write(peaks)
    else:
        processed_scan = ["{:.5f} {:.1f}".format(line[0], line[1]) for line in peaks]
        writer.write('\n'.join(processed_scan))
        writer.write("\n")

    writer.write("END IONS\n\n")
