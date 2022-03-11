import heapq
import logging
import tqdm

# Return the processed peaks as str
def padding_peaks(all_peaks, selected_peaks, thres):
    # Return experimental peaks if they are too little
    if len(all_peaks) <= thres:
        processed_scan = []
        for index in range(len(all_peaks)):
            processed_scan.append("{:.5f} {:.1f}".format(all_peaks[index][0], all_peaks[index][1]))
        return '\n'.join(processed_scan)
    
    # Pad with unlabeled peaks
    candidate_peaks = all_peaks.copy()
    assert id(candidate_peaks) != id(all_peaks)
    aux_set = set([item[0] for item in selected_peaks])
    for line in candidate_peaks:
        if line[0] in aux_set:
            candidate_peaks.remove(line)
    padding_peaks = heapq.nlargest(thres - len(selected_peaks), candidate_peaks, key=lambda x:x[1])
    assert len(selected_peaks) + len(padding_peaks) == thres, f"{len(selected_peaks)}, {len(padding_peaks)}"
    padding_peaks.extend(selected_peaks)
    padding_peaks.sort(key=lambda x:x[0])

    # Formalize
    processed_scan = []
    for line in padding_peaks:
        processed_scan.append("{:.5f} {:.1f}".format(line[0], line[1]))
    return '\n'.join(processed_scan)


def write_mgf(path, mp_spec, filtered_peaks_dict, padding_thres):
    f = open(path, "w")
    
    logging.info("writing mgf...")
    for title, filtered_peaks in tqdm.tqdm(filtered_peaks_dict.items()):
        f.write("BEGIN IONS\n")
        f.write(f"TITLE={title}\n")
        f.write(f"CHARGE={mp_spec[title][0][1]}\n")
        f.write(f"RTINSECONDS={mp_spec[title][0][2]}\n")
        f.write(f"PEPMASS={mp_spec[title][0][3]}\n")
        f.write(padding_peaks(mp_spec[title][1], filtered_peaks, padding_thres))
        f.write("\nEND IONS\n")
    f.close()    
