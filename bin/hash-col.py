#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Appends a column to a file with the hash of the selected columns in the file

Examples
--------
Examples usage::

    # hash entire lines in file
    ./hash-col.py -i NC-HAPMAP.updated.tsv | head

    # hash selected columns only under new column 'Hash'
    ./hash-col.py -i NC-HAPMAP.updated.tsv --header Hash -k CHROM POS REF ALT | head

"""
import hashlib
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

def md5_str(item):
    """
    Gets the md5sum on the string representation of an object

    Parameters
    ----------
    item:
        Python object to get the md5 sum from; should be coercible to 'str'

    Returns
    -------
    str
        the md5 hash for the item
    """
    try:
        # python 2.x
        md5 = hashlib.md5(str(item)).hexdigest()
    except:
        # python 3.x
        md5 = hashlib.md5(str(item).encode('utf-8')).hexdigest()
    return(md5)

def hash_dict(d, keys_to_hash, hash_keyname):
    """
    Creates a key in the dictionary with the hash of the given keys.

    Parameters
    ----------
    d: dict
        a dictionary to add a hash value to
    keys_to_hash: list
        a list of keys in the dictionary that should be concatenated to create a hash value. If ``None`` or empty list, all keys are used
    hash_keyname: str
        the key to add the hash value to in the dictionary

    Returns
    -------
    dict
        the original dictionary modified to include an extra key with a hash in it
    """
    if not keys_to_hash:
        hash_str = ''.join(d.values())
    else:
        hash_vals = [d[k] for k in keys_to_hash]
        hash_str = ''.join(hash_vals)
    hash_hash = md5_str(hash_str)
    d[hash_keyname] = hash_hash
    return(d)

def hash_table(fin, fout, hash_keyname, keys_to_hash, delimiter):
    """
    Hashes a tabular file, assumed to be .tsv, .csv, or similar, with headers
    """
    reader = csv.DictReader(fin, delimiter = delimiter)
    output_fieldnames = [name for name in reader.fieldnames]
    output_fieldnames.append(hash_keyname)
    writer = csv.DictWriter(fout, delimiter = delimiter, fieldnames = output_fieldnames)
    writer.writeheader()
    for row in reader:
        row = hash_dict(d = row, keys_to_hash = keys_to_hash, hash_keyname = hash_keyname)
        writer.writerow(row)

def hash_lines(fin, fout, delimiter):
    """
    Hashes the contents of a file without headers
    """
    for line in fin:
        parts = line.strip().split(delimiter)
        hash_str = ''.join(parts)
        hash_hash = md5_str(hash_str)
        fout.write(line + hash_hash + '\n')

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    delimiter = kwargs.pop('delimiter', '\t')
    header = kwargs.pop('header', None)
    keys_to_hash = kwargs.pop('keys_to_hash')


    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if header:
        hash_keyname = header
        hash_table(fin, fout, hash_keyname, keys_to_hash, delimiter)
    else:
        hash_lines(fin, fout, delimiter)

    fout.close()
    fin.close()

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-d", default = '\t', dest = 'delimiter', help="Delimiter")
    parser.add_argument("--header", default = None, dest = 'header', help="Header for the new column")
    parser.add_argument("-k", "--keys", nargs='*', dest = 'keys_to_hash', help="Column names to be hashed")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
