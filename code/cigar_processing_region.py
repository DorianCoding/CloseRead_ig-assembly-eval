#!/usr/bin/env python3

import pysam
import sys
import numpy as np
from intervaltree import Interval, IntervalTree
import argparse
import os
import warnings

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
	From https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
def calculate_mismatches(read):
    """ Calculate the number of mismatches using NM tag and CIGAR string. """
    cigartest = read.get_cigar_stats()
    if len(cigartest) == 0:
        print(f"No read for {read}, skipped")
        return -1, 0, 0, 0, 0
    cigartest= cigartest[0]
    cigarresults=read.cigartuples
    if not cigarresults:
        print(f"No read for {read}, skipped")
        return -1, 0, 0, 0, 0
    """
    nm_tag = read.get_tag('NM') if read.has_tag('NM') else 0
    insertions = sum(length for op, length in read.cigartuples if op == 1)  # Sum of 'I'
    deletions = sum(length for op, length in read.cigartuples if op == 2)  # Sum of 'D'
    ambiguous = sum(length for op, length in read.cigartuples if op == 3)  # Sum of 'N'
    soft_clipping = sum(length for op, length in read.cigartuples if op == 4)  # Sum of 'S'
    hard_clipping = sum(length for op, length in read.cigartuples if op == 5)  # Sum of 'H'
    """
    nm_tag = cigartest[-1]
    insertions = cigartest[1]  # Sum of 'I'
    deletions = cigartest[2]  # Sum of 'D'
    ambiguous = cigartest[3]  # Sum of 'N'
    soft_clipping = cigartest[4]  # Sum of 'S'
    hard_clipping = cigartest[5]  # Sum of 'H'
    mismatches = nm_tag - insertions - deletions - ambiguous
    longindels = 0
    total_indel_length = 0
    for operation, length in cigarresults:
        # Insertion or Deletion longer than 2
        if (operation == 1 or operation == 2) and length > 2:
            longindels += 1
            total_indel_length += length
    return mismatches, longindels, total_indel_length, soft_clipping, hard_clipping

def process_bam_file(bam_file_path, regions, output_dir): #, region_list_TRA, region_list_TRB, region_list_TRG):
    """ Process a BAM file to estimate mismatches for each read. """
    # Output file names
    if not os.path.exists(output_dir):
        print(f"{output_dir} was created.")
        os.makedirs(output_dir)
    else:
        print(f"{output_dir} already exists and will be overwritten.")
    output_file_names = [os.path.join(output_dir, f"{key}.txt") for key in regions.keys()]#, os.path.join(output_dir, "TRA.txt"), os.path.join(output_dir, "TRB.txt"), os.path.join(output_dir, "TRG.txt")]
    #non_overlap_file_name = os.path.join(output_dir,"nonIG.txt")

    # Open output files
    output_files = [open(name, "w") for name in output_file_names]
    #non_overlap_file = open(non_overlap_file_name, "w")

    trees = {}
    for locus, regionbig in regions.items():
        for region in regionbig:
            # Split the string into git hromosome and the 'start-end' part
            chrom, positions = region.split(':')
            # Further split the 'start-end' part into start and end positions
            start, end = positions.split('-')
            trees[locus] = [(chrom,IntervalTree())]
            if chrom != '':
               trees[locus][-1][1].addi(int(start), int(end))
            else:
               trees[locus][-1][1].addi(0, 1)
    print(trees)
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    if not bamfile.check_index():
        print("No index found, exiting")
        return
    print("")
    for i, (locus,tree_dict) in enumerate(trees.items()):
        output_files[i].write("#Read_name\tChromosome\tStart\tRead_length\tMapQ\tMismatches\tMismatch_rate\tLongindels\tTotal indel length\tIndel rate\tSoft clipping\tHard clipping\n")
        for (chr,pos) in tree_dict:
            reads = bamfile.fetch(None,None,None,str(chr) + ":" + str(pos.begin()) + "-" + str(pos.end()))
            if not reads:
                print(f"No reads found in region: {chr}:{pos}, skipping")
                continue
            for read in reads:
                print(read.get_aligned_pairs(False,True))
                print(read.get_overlap(pos.begin(),pos.end()))
                sys.exit(1)
                mismatches, longindels, total_indel_length, soft_clipping, hard_clipping = calculate_mismatches(read)
                if mismatches == -1:
                    continue #Cancel no cigar
                read_length = read.query_length if read.query_length else read.infer_query_length()
                if not read_length:
                    print(f"Read {read_name} is skipped because no length found.")
                    continue #Cannot get length so skip
                mismatch_rate = mismatches / max(1,read_length) if mismatches != 0 else 0
                read_name = read.query_name
                chromosome = bamfile.get_reference_name(read.reference_id)
                start = read.reference_start
                end = read.reference_end
                mapping_quality = read.mapping_quality
                indel_rate = total_indel_length / max(1,read_length)
                output_files[i].write(f"{read_name}\t{chromosome}\t{start}\t{read_length}\t{mapping_quality}\t{mismatches}\t{mismatch_rate}\t{longindels}\t{total_indel_length}\t{indel_rate}\t{soft_clipping}\t{hard_clipping}\n")
        printProgressBar (i+1,len(trees),"BAM analysis: ")
    '''
    for read in bamfile:
        if not read.is_unmapped and read.query_length>0:
            i+=1
            mismatches, longindels, total_indel_length, soft_clipping, hard_clipping = calculate_mismatches(read)
            read_length = read.query_length
            printProgressBar(i, bamfile.mapped, prefix = 'BAM analysis: ')
            if mismatches != 0:
                mismatch_rate = mismatches / read_length
            else:
                mismatch_rate = 0
            read_name = read.query_name
            chromosome = bamfile.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            mapping_quality = read.mapping_quality
            indel_rate = total_indel_length / read_length
            # Check for overlap with each BED file's intervals
            found_overlap = False  # Flag to check if an overlap was found
            for i, tree_dict in enumerate(trees):
                if chromosome in tree_dict and tree_dict[chromosome].overlaps(start, end):
                    output_files[i].write(f"{read_name}\t{chromosome}\t{start}\t{read_length}\t{mapping_quality}\t{mismatches}\t{mismatch_rate}\t{longindels}\t{total_indel_length}\t{indel_rate}\t{soft_clipping}\t{hard_clipping}\n")
                    found_overlap = True
                    break
            # If no overlap was found, write to the non-overlapping file
            if not found_overlap:
                non_overlap_file.write(f"{read_name}\t{chromosome}\t{start}\t{read_length}\t{mapping_quality}\t{mismatches}\t{mismatch_rate}\t{longindels}\t{total_indel_length}\t{indel_rate}\t{soft_clipping}\t{hard_clipping}\n")
    '''
    bamfile.close()
    for file in output_files:
        file.close()
    #non_overlap_file.close()


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process Cigar String..')

    # Required arguments
    parser.add_argument('input_file', help='Input SAM or BAM file.')
    parser.add_argument('IG_region', help='IG position file')
    parser.add_argument('species', help='Species name.')
    parser.add_argument('output', help='Output directory path.')

    # Parse arguments
    args = parser.parse_args()

    # Validate the input_file extension
    valid_extensions = ['.sam', '.bam']
    file_extension = os.path.splitext(args.input_file)[1]
    if file_extension not in valid_extensions:
        parser.error("Input file must be a SAM or BAM file.")

    # Initialize the lists for each gene type
    regions = {}

    # Open the file and read line by line
    with open(args.IG_region, 'r') as file:
        for line in file:
            # Split each line into its components
            parts = line.split()
            # Extract relevant data
            altbool = True if parts[1]=="alternate" else "primary"
            gene_type = parts[2]
            chr_name = parts[3]
            start = parts[4]
            end = parts[5]
            region = f"{chr_name}:{start}-{end}"
            # Append the region to the corresponding list based on gene type
            gene_type=gene_type.upper()
            if gene_type not in ["IGH","IGK","IGL","TRA","TRD","TRG","TRB"]:
                print(f"{gene_type} is not recognized as valid and is skipped.")
                continue
            elif gene_type in regions: #Haplotype
                elem=len(regions[gene_type])
                if elem>=2:
                    print(f"More than 2 regions for the locus {gene_type}, skipped.")
                    continue
                else:
                    elem.append(regions)
            else:
                regions[gene_type]=[region]


    # Output the lists to check
    print("Regions:", regions)
    outdir=os.path.join(args.output,args.species)
    if args.input_file.endswith('.sam') or args.input_file.endswith('.bam'):
        process_bam_file(args.input_file, regions, outdir) #region_list_TRA, region_list_TRB, region_list_TRG)
    else:
        raise ValueError("Not implemented.")

if __name__ == "__main__":
    main()
    
