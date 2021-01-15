
def Count_Codons(sequence):
    """takes a genetic sequence and returns a dict with the count for each of the codons"""

    __CODONS_DICT__ = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
        'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
        'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
        'TAA': 0, 'TAG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
        'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
        'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ACT': 0, 'ACC': 0,
        'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
        'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
        'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
    
    sequence = sequence.upper()
    codonsDict = __CODONS_DICT__.copy()
    for i in range(0, len(sequence), 3):
        c = sequence[i:i + 3]
        if c in codonsDict:
            codonsDict[c] += 1
    return codonsDict


def RSCU_calc(codonsDict):
    """takes a dict with the count for each of the codons present in the sequence
    and returns a dict with the RSCU values for each codon"""

    __SynonymousCodons_DICT__ = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC']}

    rscu_dict = {}
    for aa in __SynonymousCodons_DICT__:
        total = 0.0
        rscu = []
        codons = __SynonymousCodons_DICT__[aa]
        for c in codons:
            total += codonsDict[c]
            if total != 0:
                for c in codons:
                    denom = float(total) / len(codons)
                    rscu = round(codonsDict[c] / denom, 3)
                    rscu_dict[c] = rscu
    return (rscu_dict)


def Nc_calculator(seq):

    """takes a genetic sequence and returns the Nc"""

    __CODONS_DICT__ = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
        'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
        'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
        'TAA': 0, 'TAG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
        'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
        'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ACT': 0, 'ACC': 0,
        'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
        'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
        'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    __SynonymousCodons_DICT__ = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC']}


    seq=str(seq.upper())

    codonsDict=__CODONS_DICT__.copy()
    for i in range(0, len(seq), 3):
        c = seq[i:i + 3]
        if c in codonsDict:
           codonsDict[c] += 1
        else:
          raise TypeError("illegal codon:",c,"in gene")

    NcList=[]
    NcDict={}
    for aa in __SynonymousCodons_DICT__:
        P = 0.0
        n = 0.0
        codons = __SynonymousCodons_DICT__[aa]
        for c in codons:
            n += codonsDict[c]
        try:
            for c in codons:
                ni = codonsDict[c]
                Pi = float(ni / n)
                Pi = round(Pi, 3)
                P = P + (Pi * Pi)
            P = round(P, 3)
            Faa = ((n * P) - 1) / (n - 1)
            Faa = (round(Faa, 3))
            Nc = (1 / Faa)
            Nc = (round(Nc, 3))

        except ZeroDivisionError:
            Nc = 0.0
        NcList.append(Nc)
        NcDict[aa] = Nc

    return sum(NcList)


def Ncaa_calc_dict(seq):

    """takes a genetic sequence and returns a dict with the Nc for each aminoacid"""
    
    __CODONS_DICT__ = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
        'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
        'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
        'TAA': 0, 'TAG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
        'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
        'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ACT': 0, 'ACC': 0,
        'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
        'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
        'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    __SynonymousCodons_DICT__ = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC']}


    seq=str(seq.upper())

    codonsDict=__CODONS_DICT__.copy()
    for i in range(0, len(seq), 3):
        c = seq[i:i + 3]
        if c in codonsDict:
           codonsDict[c] += 1
        else:
          raise TypeError("illegal codon:",c,"in gene")

    NcList=[]
    NcDict={}
    for aa in __SynonymousCodons_DICT__:
        P = 0.0
        n = 0.0
        codons = __SynonymousCodons_DICT__[aa]
        for c in codons:
            n += codonsDict[c]
        try:
            for c in codons:
                ni = codonsDict[c]
                Pi = float(ni / n)
                Pi = round(Pi, 3)
                P = P + (Pi * Pi)
            P = round(P, 3)
            Faa = ((n * P) - 1) / (n - 1)
            Faa = (round(Faa, 3))
            Nc = (1 / Faa)
            Nc = (round(Nc, 3))

        except ZeroDivisionError:
            Nc = 0.0
        NcList.append(Nc)
        NcDict[aa] = Nc

    return NcDict

def perc_of_ala(aa_seq):

    """calculates the percentage content of alanine aa in a protein sequence"""

    aa_seq = str(aa_seq.upper())
    length = len(aa_seq)
    ala_count = 0
    for aa in aa_seq:
        if aa == "A":
            ala_count += 1
    return round((ala_count/length)*100,2)

def perc_of_cys(aa_seq):

    """calculates the percentage content of cisteine aa in a protein sequence"""
    
    aa_seq = str(aa_seq.upper())
    length = len(aa_seq)
    cys_count = 0
    for aa in aa_seq:
        if aa == "C":
            cys_count += 1
    return round((cys_count/length)*100,2)



def mannwhitney(df, variable, verbose=False):

    # Adapted from https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
    
    import pandas as pd
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # seed the random number generator
    seed(1)

    # low and high
    selection = [variable, 'Expression_Level']
    df2 = df[selection]
    high = df[df.Expression_Level == 'high']
    high = high[variable]

    selection = [variable, 'Expression_Level']
    df2 = df[selection]
    low = df[df.Expression_Level == 'low']
    low = low[variable]

    # compare samples
    stat, p = mannwhitneyu(low, high)
    # print('Statistics=%.3f, p=%.3f' % (stat, p))

    # interpret
    alpha = 0.05

    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'

    results = pd.DataFrame({'variable': variable,
                            'Statistics': stat,
                            'p': p,
                            'alpha': alpha,
                            'Interpretation': interpretation}, index=[0])
    filename = 'mannwhitneyu_' + variable + "_" + '.csv'
    results.to_csv("statistical_tests/{}".format(filename))

    return results

if __name__ == "__main__":
    Count_Codons()
    RSCU_calc()
    displayOutput()
    Nc_calculator()
    perc_of_ala()
    perc_of_cys()
    pass