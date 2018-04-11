SHELL:=/bin/bash
REFDIR:=/ifs/data/sequence/results/external/NYU/snuderllab/ref
ANNOVAR_DB_DIR:=/ifs/data/molecpathlab/bin/annovar/db/hg19
.PHONY: annovar_db ref

# ~~~~~ SETUP PIPELINE ~~~~~ #
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

update: ./nextflow
	./nextflow self-update

annovar_db: install
	[ ! -d "$(ANNOVAR_DB_DIR)" ] && echo "> system ANNOVAR db dir does not exist, setting up local dir..." && { \
	./nextflow run annovar_db.nf -profile ref ; } || :


ref: install
	[ ! -d "$(REFDIR)" ] && echo "> system ref dir doesnt exist, setting up local ref dir..." && { \
	./nextflow run ref.nf -profile ref ; } || :

setup: annovar_db ref

# ~~~~~ RUN PIPELINE ~~~~~ #
test: install
	./nextflow run test.nf

run: install
	./nextflow run main.nf -with-dag flowchart.dot $(EP)

flowchart:
	[ -f flowchart.dot ] && dot flowchart.dot -Tpng -o flowchart.png

# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace*.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-reports:
	rm -f *.html.*

clean-flowcharts:
	rm -f *.dot.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

# deletes files from previous runs of the pipeline, keeps current results
clean: clean-logs clean-traces clean-reports clean-flowcharts

# deletes all pipeline output
clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
