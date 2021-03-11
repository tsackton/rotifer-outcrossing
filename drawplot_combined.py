import matplotlib.pyplot as plt
import re
import sys
import os
import argparse
from Bio import SeqIO
import itertools

def create_plot(filename, filter):

    #Get region name
    region_num = re.search("(\d+).*", os.path.basename(filename)).group(1)
    region_name = "Region: " + region_num

    #check to see if this is on the acceptable list:
    try:
        if (region_num not in filter):
            print("Skipping region number " + region_num)
            return;
    except:
        pass

    print("Plotting region number " + region_num)

    # Get fasta file
    seq_dict = SeqIO.index(filename, "fasta")
    seq_names = list(seq_dict.keys())

    #number of sequences
    seq_num = range(1,len(seq_dict)+1)
    pairwise = list(itertools.combinations(seq_num, 2))
    combo_num = len(pairwise)

    # Get sequences
    fig, ax = plt.subplots(combo_num)
    i = -1
    figname = re.sub("\.fa\w+", '', os.path.basename(filename)) + ".pdf"

    for seq_pair in pairwise:
        # Iterate the counter
        i = i+1
        # Select the next axis in the figure
        axis = ax[i]
        seq1_name = seq_names[seq_pair[0]-1]
        seq2_name = seq_names[seq_pair[1]-1]
        seq1 = seq_dict[seq1_name]
        seq2 = seq_dict[seq2_name]
        if (len(seq1) != len(seq2)):  # in case the sequences didn't past correctly
            print("sequences are not aligned properly")
            exit()
        mismatchL = []  # list to hold all the matches/mismatches as 0/1
        for index in range(len(seq1)):
            # assuming gaps are represented by dashes
            if (seq1[index] != '-') and (seq2[index] != '-'):
                if (seq1[index] == seq2[index]):
                    mismatchL.append(0)
                else:
                    mismatchL.append(1)

        # mismatch number
        num_mismatch = sum(mismatchL)

        # make plot
        # set plot range
        axis.set_xlim(0, len(seq1) + 1)
        axis.set_ylim(0, 1)
        axis.set_yticks([])

        # draws verticle lines of length 5 whenever there's a mismatch
        for index in range(len(mismatchL)):
            if (mismatchL[index] == 1):
                axis.vlines(index+1, 0, 1)

        axis.set_xlabel('Position')  # x axis title
        axis.set_title(seq1_name + " vs " + seq2_name +
            "\nNumber of differences: " + str(num_mismatch))  # chart title

    fig.set_size_inches(10, 10)
    fig.suptitle(region_name, fontsize=16)

    #clean up layout
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(figname)
    plt.close()

def parse_filter(filterfile):

    if (filterfile == None):
        return(None)
    else:
        filterlist = [line.rstrip('\n') for line in open(filterfile)]

    return(filterlist)

def main():
    # Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta", help="Fasta file of alignment to process, or directory of fasta files", required=True)
    parser.add_argument(
        "--filter", help="Text file of regions to keep", required=False, default=None)
    args = parser.parse_args()

    filter_list = parse_filter(args.filter)

    if (os.path.isfile(args.fasta)):
        create_plot(args.fasta, filter_list)
    elif (os.path.isdir(args.fasta)):
        all_files = [f for f in os.listdir(args.fasta) if f.endswith(('.fasta', '.fas', '.fa'))]
        for seqfile in all_files:
            seqfile = args.fasta + "/" + seqfile
            create_plot(seqfile, filter_list)
    else:
        print("You shouldn't be here")



if __name__ == "__main__":
    main()
