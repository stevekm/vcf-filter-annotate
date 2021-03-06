#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrapper script to import a table to a SQLite file

Examples
--------
Example usage::

    table2sqlite.py input_file.tsv output_file.sqlite table_name <delimiter>

    ./table2sqlite.py ../output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.updated.tsv output.sqlite variants $'\t'

Notes
-----
Does not create primary keys. Calling repeatedly will import duplicate columns. Only use this for create a new SQLite database, not for updating a previous one. Imports a single file into a single SQLite table with the given name.
"""
import sys
import sqlite3
import argparse
from util import sqlite_tools as sqt

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file')
    output_file = kwargs.pop('output_file')
    table_name = kwargs.pop('table_name')
    delimiter = kwargs.pop('delimiter', '\t')
    dump_csv = kwargs.pop('dump_csv', None)
    dump_sqlite = kwargs.pop('dump_sqlite', None)

    conn = sqlite3.connect(output_file)
    sqt.csv2sqlite(conn = conn, input_file = input_file, table_name = table_name, delimiter = delimiter)

    if dump_csv:
        sqt.dump_csv(conn = conn, table_name = table_name, output_file = dump_csv)
    if dump_sqlite:
        sqt.dump_sqlite(conn = conn, output_file = dump_sqlite)

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", required = True, dest = 'input_file', help="Input file")
    parser.add_argument("-o", required = True, dest = 'output_file', help="Output file")
    parser.add_argument("-t", "--table-name", required = True, dest = 'table_name', help="Name for the SQLite table in the database")
    parser.add_argument("-d", "--delimiter", default = '\t', dest = 'delimiter', help="Delimiter")
    parser.add_argument("--dump-csv", default = None, dest = 'dump_csv', help="File to dump .csv formatted output to")
    parser.add_argument("--dump-sqlite", default = None, dest = 'dump_sqlite', help="File to dump SQLite database contents to")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
