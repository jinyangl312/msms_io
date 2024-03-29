import re
from msms_io.mass.ion_calc import calc_pepmass


def sort_alpha_beta_order(mixed_peptide, mixed_modification):
    '''
    Rearrange the order of alpha and beta peptide. Rearrange modification at the same time.
    Used for the comparison between search engines.
    '''

    line = re.split("\-|\(|\)", mixed_peptide)
    if len(line) != 6:
        return mixed_peptide, mixed_modification

    seq1, site1, seq2, site2 = line[0], line[1], line[3], line[4]

    def compare_mixed_modification(mixed_modification, seq1_length):
        '''
        Define a comparison function for modifications.
        Return True if the mod str of beta peptide is greater.
        '''

        if mixed_modification == "" or mixed_modification == "null":
            return False

        # Split mods
        seq1_mod_list = []
        seq2_mod_list = []
        for mod in re.split("\;", mixed_modification):
            mod_name = re.search(".+(?=\(\d+)", mod).group()
            mod_site = int(re.search("(?<=\()\d+(?=\))", mod).group())
            if mod_site <= seq1_length+1:
                seq1_mod_list.append(f"{mod_name}({mod_site})")
            else:
                seq2_mod_list.append(f"{mod_name}({mod_site-seq1_length-3})")
        # Compare mods
        return ";".join(seq1_mod_list) < ";".join(seq2_mod_list)

    # alpha peptide: has greater pepmass
    # If has the same pepmass with beta peptide, choose the one with greater link site
    # If has the same link site, compare modification string
    if calc_pepmass(seq1) < calc_pepmass(seq2) or (  # has greater pepmass
        calc_pepmass(seq1) == calc_pepmass(seq2) and (
            int(site1) < int(site2) or (  # equal pepmass but with greater link site
                int(site1) == int(site2) and (
                    compare_mixed_modification(mixed_modification, len(seq1)))))):  # but with greater modification
        # Rearrange modification for alpha and beta
        if not (mixed_modification == "" or mixed_modification == "null"):
            new_mod_list_1 = []
            new_mod_list_2 = []
            for mod in re.split("\;", mixed_modification):
                mod_name = re.search(".+(?=\(\d+)", mod).group()
                mod_site = int(re.search("(?<=\()\d+(?=\))", mod).group())
                if mod_site <= len(seq1):
                    new_mod_list_1.append(
                        f"{mod_name}({mod_site+len(seq2)+3})")
                else:
                    new_mod_list_2.append(
                        f"{mod_name}({mod_site-len(seq1)-3})")
            new_mod_list = new_mod_list_2 + new_mod_list_1
            mixed_modification = ";".join(new_mod_list)

        # Rearrange sequence and site for alpha and beta
        mixed_peptide = f'{seq2}({site2})-{seq1}({site1})'

    return mixed_peptide, mixed_modification


def assemble_alpha_beta_peptides(seq1, seq2, site1, site2, mods1, mods2):
    '''
    Assemble two peptides into a cross-linked peptide.
    Used in pFindLink
    '''

    mod_list_1 = []
    if not mods1 == "":
        for mod in re.split("\;", mods1.strip(";")):
            # 24,Deamidated[N];33,Carbamidomethyl[C];
            mod_site = int(re.search("\d+(?=\,)", mod).group())
            mod_name = re.search("(?<=\,).+(?=$)", mod).group()
            mod_list_1.append([mod_site, mod_name])
    mod_list_2 = []
    if not mods2 == "":
        for mod in re.split("\;", mods2.strip(";")):
            # 24,Deamidated[N];33,Carbamidomethyl[C];
            mod_site = int(re.search("\d+(?=\,)", mod).group())
            mod_name = re.search("(?<=\,).+(?=$)", mod).group()
            mod_list_2.append([mod_site, mod_name])

    new_mod_list = [f"{x[1]}({x[0]})" for x in mod_list_1] + \
        [f"{x[1]}({x[0]+len(seq1)+3})" for x in mod_list_2]

    return f'{seq1}({site1})-{seq2}({site2})', ";".join(new_mod_list)


def sort_site_order(mixed_line):
    '''
    Rearrange the order of cross-linked sites.
    '''

    if mixed_line == '':
        return ''

    res = []
    for mixed_site in mixed_line.strip("/").split("/"):
        pre_line = re.split('(?<=\))\-', mixed_line)
        protein_1, site1 = re.split("\(|\)", pre_line[0])[0:2]
        protein_2, site2 = re.split("\(|\)", pre_line[1])[0:2]

        if protein_1 > protein_2 or (  # string order
            protein_1 == protein_2 and (
                int(site1) > int(site2))):

            # Rearrange protein and site
            res.append(f'{protein_2}({site2})-{protein_1}({site1})/')
        else:
            res.append(f'{mixed_site}/')
    if len(res) > 1:
        res.sort()

    return "".join(res)


def sort_protein_order(mixed_line):
    '''
    Rearrange the order of cross-linked proteins.
    '''

    res = set()
    for mixed_site in mixed_line.strip("/").split("/"):
        line = re.split("\-", mixed_site)
        protein_1, protein_2 = line[0], line[1]

        if protein_1 > protein_2:  # string order
            # Rearrange protein and site
            res.add(f'{protein_2}-{protein_1}/')
        else:
            res.add(f'{mixed_site}/')

    res = list(res)
    if len(res) > 1:
        res.sort()

    return "".join(res)


def split_xl_protein(proteins_line):
    return re.sub("\-", "/", proteins_line)


def rmv_xl_site(proteins_line):
    return re.sub("\(\d+\)", "", proteins_line)


def is_intra_protein(mixed_protein):
    line = re.split("\-\/", mixed_protein)
    return line[0] == line[1]
