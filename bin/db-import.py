#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import csv
import sqlite3
from util import sqlite_tools as sqt
recalc_tsv_file = "/Users/steve/projects/vcf-filter-annotate/output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.recalc.tsv"
recalc_avinput_file = "/Users/steve/projects/vcf-filter-annotate/output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.avinput.recalc.tsv"
annotations_file = "/Users/steve/projects/vcf-filter-annotate/output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.hg19_multianno.txt"
sampleID = "NC-HAPMAP"
sequencing_run = "12345"

#  setup db
db_file = "annotations.sqlite"
conn = sqlite3.connect(db_file)
sqt.create_table(conn = conn, table_name = "annotations", col_name = "hash", col_type = "TEXT", is_primary_key = True)

# read in avinput file
with open(recalc_avinput_file) as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
        row['hash'] = sqt.md5_str(''.join([ row["Chr"], row["Start"], row["Ref"], row["Alt"], row["CHROM"], row["POS"], row["REF"], row["ALT"], sampleID, sequencing_run ]))
        row = sqt.sanitize_dict_keys(d = row)
        for key in row.keys():
            sqt.add_column(conn = conn, table_name = "annotations", col_name = key, col_type = "TEXT")
        sqt.sqlite_insert(conn = conn, table_name = "annotations", row = row)
        print(row)



# ~~~~~ CLEAN UP ~~~~~ #
# dump the entire database
db_dump_file = "annotations.sqlite.dump.txt"
print("dumping to file: {0}".format(db_dump_file))
sqt.dump_sqlite(conn = conn, output_file = db_dump_file)

# create csv dumps of database
table_names = sqt.get_table_names(conn = conn)
for name in table_names:
        output_file = "{0}.{1}.csv".format("annotations", name)
        print("dumping to file: {0}".format(output_file))
        sqt.dump_csv(conn = conn, table_name = name, output_file = output_file)
conn.commit()
conn.close()
