#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recalculates the variant allele frequency for GATK Haplotype Caller output

INPUT: table-formatted .vcf output from GATK VariantsToTable
OUTPUT: .TSV formatted table with recalculated variant allele frequency in the FREQ column
USAGE: ./recalc-vcf-AF.py NC-HAPMAP.vcf.txt NC-HAPMAP GATKHC
"""
import csv
import sys
import argparse

def GATKHC(fin, fout, sampleID):
    """
    Recalculates the variant allele frequency for GATK Haplotype Caller output
    assumes VCFv4.2

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file
    sampleID: str
        identifier for the sample in the input file connection

    """
    # column names for the AD and DP fields in the table
    AD_key = "{0}.AD".format(sampleID)
    DP_key = "{0}.DP".format(sampleID)

    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        ref_AD = float(row[AD_key].split(',')[0])
        alt_AD = float(row[AD_key].split(',')[1])
        depth = float(row[DP_key])
        ref_DP = round(ref_AD / depth, 4)
        alt_AF = round(alt_AD / depth, 6)
        row['FREQ'] = alt_AF
        writer.writerow(row)

def LoFreq(fin, fout):
    """
    LoFreq does not need recalulating; output a new column in the file with 'FREQ' for consistency

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    # allele frequency column
    AF_key = "AF"
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        AF_val = row[AF_key]
        row['FREQ'] = AF_val
        writer.writerow(row)

def MuTect2(fin, fout):
    """
    Outputs the TUMOR.AF as a new column called 'FREQ', and outputs the TLOD value as QUAL

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        # set the FREQ to the tumor AF
        tumor_AF_value = round(row[tumor_AF_key], 6)
        row['FREQ'] = row['TUMOR.AF']
        # change QUAL to the TLOD
        row['QUAL'] = row['TLOD']
        writer.writerow(row)


def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    caller = kwargs.pop('caller')
    sampleID = kwargs.pop('sampleID')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if caller == "GATKHC":
        GATKHC(fin, fout, sampleID)
        fout.close()
        fin.close()
    elif caller == "LoFreq":
        LoFreq(fin, fout)
        fout.close()
        fin.close()
    elif caller == "MuTect2":
        MuTect2(fin, fout)
        fout.close()
        fin.close()
    else:
        print("ERROR: caller not recognized: {0}".format(caller))
        sys.exit(1)




def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")

    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    parser.add_argument("-s", "--sampleID", dest = 'sampleID', help="Sample ID", required=True)
    args = parser.parse_args()

    main(**vars(args))



if __name__ == '__main__':
    parse()
