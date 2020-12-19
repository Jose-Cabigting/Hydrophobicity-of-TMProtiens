import matplotlib.pyplot as plt
from argparse import ArgumentParser

#some of the functions have similar lines of code with small changes to show different outputs I'm still trying to simplify

def seq_score_conv(seq, hydro_dictionary):
    """ Inputs are Amino acid sequence in all Caps and a hydrophobicity scale dictionary,
    Outputs are a table of hydrophobicity scores """
    hydo_score_conv = []
    for aa in seq:
        if (aa in hydro_dictionary):
            hydo_score_conv.append(hydro_dictionary[aa])
    return hydo_score_conv

def hydro_sum(seq, hydro_dictionary):
    """ Inputs are Amino acid sequence and hydrophobicity dictionary,
    Outputs are a sum of the numerical values from the list of scores created """
    hydroNum = 0
    for aa in seq:
        if aa in hydro_dictionary:
            num = float(hydro_dictionary[aa])
            hydroNum += num
            hydroNum = hydroNum
    return hydroNum

def stepbystep_hydro_sum(seq, hydro_dictionary):
    """ Inputs are Amino acid sequence and hydrophobicity dictionary,
    Output shows a list of numerical values in which each score is added step by step """
    sum_steps = []
    hydroNum = 0
    for aa in seq:
        if aa in hydro_dictionary:
            num = float(hydro_dictionary[aa])
            hydroNum += num
            sum_steps.append(hydroNum)
            hydroNum = hydroNum
    return sum_steps

def average_hydro_score(int, seq):
    """ Inputs are an integer value (the hydrophobicity score sum) and the original Amino acid sequence,
    Output is an average numerical value to the third decimal """
    sum = int
    seq_length = len(seq)
    average_score = str(round(sum/seq_length, 3))
    return average_score

def hydro_plot(list):
    """Input is a list of hydophobicity scores that are added together step by step,
    Output is a line graph showing hydrophobicity score over length of sequence"""
    plt.plot(list)
    plt.ylabel("Hydrophobicity Values")
    plt.xlabel("Amino Acid location (starting from position 0)")
    plt.show()

def hydro_viz(seq, hydro_dictionary, seperator=''):
    """Input is an amino acid sequence and a hydrophobicity dictionary along with a seperator/joining chracter,
    Output is a string of chracters (+ or -) that represenet phobic or philic tendencies"""
    list_viz = []
    for aa in seq:
        if aa in hydro_dictionary:
            num = float(hydro_dictionary[aa])
            if num > 0:
                list_viz.append('+')
            else:
                list_viz.append('-')
    return seperator.join(list_viz)

def protien_pass_count(seq):
    """ Input is a string or sequence of + or - chracters, Output is a numerical value counting "+-" in the seq """
    plus_minus_count = seq.count("+-")
    return plus_minus_count

def passage_comparison(int1, int2):
    """ Input is two numerical values a known number of times the protien passes the membrane and the found,
    Output is a true or false statement, values match or values don't match """
    found_pass_number = int(int1)
    known_pass_number = int(int2)
    match = ("Found and Known Values are Equal")
    miss = ("Found and Known Values are Not Equal")
    if found_pass_number == known_pass_number:
        return match
    else:
        return miss

kyte_doolittle_hydro_scale = {
#This Hydrophobicity scale is one of many, if one would like to change to a different scale then alter the values for the AA's
    'I':'4.5',
    'F':'2.8',
    'V':'4.2',
    'L':'3.8',
    'W':'-0.9',
    'M':'1.9',
    'A':'1.8',
    'G':'-0.4',
    'C':'2.5',
    'Y':'-1.3',
    'P':'-1.6',
    'T':'-0.7',
    'S':'-0.8',
    'H':'-3.2',
    'E':'-3.5',
    'N':'-3.5',
    'Q':'-3.5',
    'D':'-3.5',
    'K':'-3.9',
    'R':'-4.5'
    }
    
parser = ArgumentParser()
parser.add_argument("--seq", help="intended amino acid sequence submitted for analysis")
parser.add_argument("--cross", help="number of known times transmembrane protien passes through cell membrane")
args = parser.parse_args()

seq = args.seq
cross_number = args.cross

example_seq = f"{seq}"
practice_seq = example_seq.upper()
#This is where the Amino acid sequnece being analyzed witll be inserted, file reader function pending
example_known = f"{cross_number}"
practice_known = int(example_known)
#This known value of number of times the transmembrane eprotien passes the membrane is inputed here

practice_conv = seq_score_conv(practice_seq, kyte_doolittle_hydro_scale)
practice_sum = hydro_sum(practice_seq, kyte_doolittle_hydro_scale)
practice_stepbystep = stepbystep_hydro_sum(practice_seq,kyte_doolittle_hydro_scale)
practice_seq_length = len(practice_seq)
practice_avg = average_hydro_score(practice_sum, practice_seq)
# avg_hydro_score = str(round(practice_sum/practice_seq_length, 3)) alternative way to achieve average
practice_viz = hydro_viz(practice_seq,kyte_doolittle_hydro_scale)
practice_count = protien_pass_count(practice_viz)
practice_comp = passage_comparison(practice_count, practice_known)

print ("Sequence Submitted for Analysis:", practice_seq)
print ("Individual Hyrdophobicity scores:", practice_conv)
print ("Step by Step Summation:", practice_stepbystep)
print ("Total Hydrophobicity Score:", practice_sum)
print ("Average Hydrophobicity Score:", practice_avg)
# print ("Average Hydrophobicity score:", avg_hydro_score) alternative representation of calculated average
print ("Hydrophobic (+) or Hydrophilic (-) Visualization:", practice_viz)
print ("# of Times Cell Passes Through Membrane (Found):", practice_count)
print ("# of Times Cell Passes Through Membrane (Known):", practice_known)
print (practice_comp)

hydro_plot(practice_stepbystep)
