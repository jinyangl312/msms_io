import logging


def parse_peaks(all_peaks, theoretical_peaks, deviation=20e-6):
    '''
    Given two sorted peaks, parse them with deviation in O(n) complexity.
    Allow multiple theoretical peaks match one peak. Do not allow
    multiple peaks for one theoretical peak, but choose the greater peak
    if there is another peak within the window (only consider two peaks).
    Used in the build of peaks alignment table.
    '''

    selected_peaks = []
    peaks_id = 0
    candidate = all_peaks[peaks_id]
    # For every theoretical peak
    for theo_peak in theoretical_peaks:
        # Continue if the theo peak is too small than the candidate
        if (candidate[0] / theo_peak[0] - 1) > deviation:
            continue

        # update candidate if it is too small than the theo peak
        while (candidate[0] / theo_peak[0] - 1) < -deviation:
            peaks_id += 1
            # if parse completed
            if peaks_id == len(all_peaks):
                return selected_peaks
            candidate = all_peaks[peaks_id]

        # Continue if the theo peak is too small than the candidate again
        # (Or the candidate is getting too large than the theo peak)
        if (candidate[0] / theo_peak[0] - 1) > deviation:
            continue

        # Now the two peaks should be within 20 ppm
        # Append results
        if peaks_id < len(all_peaks)-1 and abs(all_peaks[peaks_id+1][0] / theo_peak[0] - 1) < deviation:
            logging.warning(
                f'Multiple peaks for one theoretical peak: {theo_peak}\n1. {candidate}\n2. {all_peaks[peaks_id+1]}\n')
            # Match the peak with greater intensity
            if all_peaks[peaks_id+1][1] > all_peaks[peaks_id][1]:
                selected_peaks.append((theo_peak[1], all_peaks[peaks_id+1][1]))
            else:
                selected_peaks.append((theo_peak[1], candidate[1]))
        else:
            selected_peaks.append((theo_peak[1], candidate[1]))

        #peaks_id += 1
        if peaks_id == len(all_peaks):
            return selected_peaks
        candidate = all_peaks[peaks_id]

    return selected_peaks


def parse_peaks_with_suffix_tag(all_peaks, tag_length, deviation=20e-6):
    '''
    Given a list of sorted peaks, find the peaks with a suffix tag
    with deviation in O(n) complexity.
    Forked from parse_peaks.
    Used in ?
    '''

    selected_peaks = []
    peaks_id = 0
    candidate2 = all_peaks[peaks_id]
    # For every theoretical peak
    for candidate1 in all_peaks:
        # Continue if the theo peak is too small than the candidate
        if (1 - (candidate1[0] + tag_length) / candidate2[0]) > deviation:
            continue

        # update candidate if it is too small than the theo peak
        while (1 - (candidate1[0] + tag_length) / candidate2[0]) < -deviation:
            peaks_id += 1
            # if parse completed
            if peaks_id == len(all_peaks):
                return selected_peaks
            candidate2 = all_peaks[peaks_id]

        # Continue if the theo peak is too small than the candidate again
        # (Or the candidate is getting too large than the theo peak)
        if (1 - (candidate1[0] + tag_length) / candidate2[0]) > deviation:
            continue

        # Now the two peaks should be within 20 ppm
        # Append results
        selected_peaks.append((candidate1, candidate2))
        if peaks_id < len(all_peaks)-1 and \
                abs(all_peaks[peaks_id+1][0] / (candidate1[0] + tag_length) - 1) < deviation:
            logging.warning(
                f'Multiple peaks for one theoretical peak: {candidate1}\n{candidate2}\n{all_peaks[peaks_id+1]}')
        # peaks_id += 1 # Multiple theo peaks for one peak
        if peaks_id == len(all_peaks):
            return selected_peaks
        candidate2 = all_peaks[peaks_id]

    return selected_peaks


def parse_peaks_for_filter(all_peaks, theoretical_peaks):
    '''
    Forked from parse_peaks.
    Used in ?
    '''
    selected_peaks = []
    deviation = 20e-6

    peaks_id = 0
    candidate = all_peaks[peaks_id]
    # For every theoretical peak
    for theo_peak in theoretical_peaks:
        # Continue if the theo peak is too small than the candidate
        if (candidate[0] / theo_peak[0] - 1) > deviation:
            continue

        # update candidate if it is too small than the theo peak
        while (candidate[0] / theo_peak[0] - 1) < -deviation:
            peaks_id += 1
            # if parse completed
            if peaks_id == len(all_peaks):
                return selected_peaks
            candidate = all_peaks[peaks_id]

        # Continue if the theo peak is too small than the candidate again
        # (Or the candidate is getting too large than the theo peak)
        if (candidate[0] / theo_peak[0] - 1) > deviation:
            continue

        # Now the two peaks should be within 20 ppm
        # Append results
        selected_peaks.append((candidate))
        if abs(all_peaks[peaks_id+1][0] / theo_peak[0] - 1) < deviation:
            logging.warning(
                f'Multiple peaks for one theo peak: {theo_peak}\n{candidate}\n{all_peaks[peaks_id]}')
        # peaks_id += 1 # Multiple theo peaks for one peak
        if peaks_id == len(all_peaks):
            return selected_peaks
        candidate = all_peaks[peaks_id]

    return selected_peaks
