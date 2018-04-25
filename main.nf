// ~~~~~~~~~~ SETUP PARAMETERS ~~~~~~~~~~ //
params.runID = "170519_NB501073_0010_AHCLLMBGX2"
params.resultsID = null

// set a timestamp variable if resultsID not passed
import java.text.SimpleDateFormat
def resultsID
if ( params.resultsID == null ) {
    Date now = new Date()
    SimpleDateFormat timestamp = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss")
    resultsID = timestamp.format(now)
} else {
    resultsID = params.resultsID
}

// path to the output directory
output_dir_path = new File(params.output_dir).getCanonicalPath()

// path to the current directory
current_dir_path = new File(System.getProperty("user.dir")).getCanonicalPath()

// get the system hostname to identify which system the pipeline is running from
String localhostname = java.net.InetAddress.getLocalHost().getHostName();

println ">>> ${params.runID} ${resultsID} ${output_dir_path} ${current_dir_path} ${localhostname}"

// ~~~~~~~~~~ INPUT CHANNELS ~~~~~~~~~~ //
Channel.from( [
            ['HaplotypeCaller', 'NC-HAPMAP', file("input/NC-HAPMAP.GATKHC.vcf.gz")],
            ['LoFreq', 'NC-HAPMAP', file("input/NC-HAPMAP.LoFreq.vcf.gz")],
            ['HaplotypeCaller', 'SC-SERACARE', file("input/SC-SERACARE.GATKHC.vcf.gz")],
            ['LoFreq', 'SC-SERACARE', file("input/SC-SERACARE.LoFreq.vcf.gz")],
            [ "LoFreq", "NTC-H2O", file("input/vcf_lofreq/NTC-H2O.vcf.gz") ],
            [ "LoFreq", "HapMap-B17-1267", file("input/vcf_lofreq/HapMap-B17-1267.vcf.gz") ],
            [ "LoFreq", "SeraCare-1to1-Positive", file("input/vcf_lofreq/SeraCare-1to1-Positive.vcf.gz") ],
            [ "HaplotypeCaller", "NTC-H2O", file("input/vcf_hc/NTC-H2O.vcf.gz") ],
            [ "HaplotypeCaller", "HapMap-B17-1267", file("input/vcf_hc/HapMap-B17-1267.vcf.gz") ],
            [ "HaplotypeCaller", "SeraCare-1to1-Positive", file("input/vcf_hc/SeraCare-1to1-Positive.vcf.gz") ]
            ] )
        .set { sample_variants_zipped }

Channel.from([
    // caller, tumorID, normalID, chrom, file
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr8", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr8.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr2", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr2.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr5", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr5.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr12", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr12.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr3", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr3.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr10", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr10.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr20", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr20.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr14", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr14.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chrX", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrX.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr13", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr13.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr19", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr19.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr22", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr22.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr21", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr21.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr15", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr15.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr9", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr9.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr4", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr4.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr7", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr7.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr16", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr16.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr17", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr17.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr18", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr18.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr11", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr11.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr1", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr1.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chrY", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrY.vcf.gz") ],
[ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "chr6", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr6.vcf.gz") ]
    ])
    .map { caller, tumorID, normalID, chrom, input_file ->
        def pairID = "${tumorID}.${normalID}.${chrom}"
        return [ caller, pairID, input_file ]
    }
    .set { pairs_variants_zipped }


Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
        .into { ref_fasta; ref_fasta2; ref_fasta3 }

Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")
        .into { ref_fai; ref_fai2; ref_fai3 }
Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict")
        .into { ref_dict; ref_dict2; ref_dict3 }
Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }


// ~~~~~~~~~~ ANALYSIS PIPELINE ~~~~~~~~~~ //
// # 1. left normalize indels & split multiallelic entries (.vcf -> norm.vcf)
// # 2. filtere vcf (norm.vcf -> .norm.filtered.vcf)
// # 2. convert to tsv (.norm.filtered.vcf -> .tsv)
// # 3. recalculate allele frequency (.tsv -> reformat.tsv)
// # 4. annotate (norm.vcf -> .avinput, .hg19_multianno.txt)
// # 5. merge annotations with .tsv (.hg19_multianno.txt, reformat.tsv -> ... )

process unzip_samples {
    // unzip the vcf.gz files
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/unzip_samples", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf) from sample_variants_zipped.concat(pairs_variants_zipped)

    output:
    set val(caller), val(sampleID), file("${prefix}.vcf") into sample_variants

    script:
    prefix = "${sampleID}.${caller}"
    """
    gunzip -c "${sample_vcf}" > "${prefix}.vcf"
    """
}

process normalize_vcf {
    // normalize and split the VCF entries
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/normalize_vcf", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta) from sample_variants.combine(ref_fasta)

    output:
    set val(caller), val(sampleID), file("${prefix}.norm.vcf") into (normalized_variants, normalized_variants2)
    file("${prefix}.bcftools.multiallelics.stats.txt")
    file("${prefix}.bcftools.realign.stats.txt")

    script:
    prefix = "${sampleID}.${caller}"
    """
    cat ${sample_vcf} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${prefix}.bcftools.multiallelics.stats.txt" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${prefix}.bcftools.realign.stats.txt" > \
    "${prefix}.norm.vcf"
    """
}

process filter_vcf {
    // filter the VCF entries
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/filter_vcf", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from normalized_variants2.combine(ref_fasta2).combine(ref_fai2).combine(ref_dict2)

    output:
    set val(caller), val(sampleID), file("${prefix}.filtered.vcf") into filtered_vcfs

    script:
    prefix = "${sampleID}.${caller}"
    if( caller == 'HaplotypeCaller' )
        """
        # report if
        # alternate allele freq (allele depth / depth) greater than 0.5
        # more than 5 variant call supporting reads
        # quality reads present (reported depth >0)
        gatk.sh -T SelectVariants \
            -R "${ref_fasta}" \
            -V "${sample_vcf}" \
            --sample_name "${sampleID}" \
            -select "vc.getGenotype('${sampleID}').getAD().1 / vc.getGenotype('${sampleID}').getDP() > 0.50" \
            -select "vc.getGenotype('${sampleID}').getAD().1 > 5" \
            -select "vc.getGenotype('${sampleID}').getDP() > 0" \
            > "${prefix}.filtered.vcf"
        """
    else if( caller == 'LoFreq' )
        """
        # do not report if frequency is less than 1%
        gatk.sh -T SelectVariants \
                -R "${ref_fasta}" \
                -V "${sample_vcf}" \
                -select "AF > 0.01"  \
                > "${prefix}.filtered.vcf"
        """
    else if( caller == 'MuTect2' )
        """
        # report if
        # T frequency is more than 3%
        # N frequency is less than 5%
        # at least 5 variant call supporting reads
        # T frequency is sufficiently higher than N frequency # "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
        # only 'PASS' entries
        gatk.sh -T SelectVariants \
                -R "${ref_fasta}" \
                -V "${sample_vcf}" \
                -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) )  > 0.03" \
                -select "(vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) )  > 0.05" \
                -select "vc.getGenotype('TUMOR').getAD().1 > 5" \
                -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) ) > (vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) ) * 5" \
                -select 'vc.isNotFiltered()' \
                > "${prefix}.filtered.vcf"
        """
    else
        error "Invalid caller: ${caller}"
}

process vcf_2_tsv {
    // convert VCF file to TSV
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/vcf_2_tsv", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcfs.combine(ref_fasta3).combine(ref_fai3).combine(ref_dict3)

    output:
    set val(caller), val(sampleID), file(sample_vcf), file("${prefix}.tsv") into vcf_tsvs

    script:
    prefix = "${sampleID}.${caller}"
    if( caller == 'HaplotypeCaller' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
        -GF AD -GF DP \
        -o "${prefix}.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F SB -F INDEL -F CONSVAR -F HRUN \
        -o "${prefix}.tsv"
        """
    else if( caller == 'MuTect2' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
        -GF AD -GF DP -GF AF \
        -o "${prefix}.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process reformat_vcf_tsv {
    // reformat and adjust the TSV table for consistency downstream
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/reformat_vcf_tsv", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(sample_tsv) from vcf_tsvs

    output:
    set val(caller), val(sampleID), file(sample_vcf), file("${prefix}.reformat.tsv") into vcfs_tsvs_reformat

    script:
    prefix = "${sampleID}.${caller}"
    if( caller == 'HaplotypeCaller' )
        """
        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c GATKHC -s "${sampleID}" -i "${sample_tsv}" | \
        paste-col.py --header "Sample" -v "${sampleID}"  | \
        paste-col.py --header "Run" -v "${params.runID}" | \
        paste-col.py --header "Results" -v "${resultsID}" | \
        paste-col.py --header "Location" -v "${current_dir_path}" | \
        paste-col.py --header "VariantCaller" -v "${caller}" | \
        paste-col.py --header "System" -v "${localhostname}" > \
        "${prefix}.reformat.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        reformat-vcf-table.py -c LoFreq -s "${sampleID}" -i "${sample_tsv}" | \
        paste-col.py --header "Sample" -v "${sampleID}"  | \
        paste-col.py --header "Run" -v "${params.runID}" | \
        paste-col.py --header "Results" -v "${resultsID}" | \
        paste-col.py --header "Location" -v "${current_dir_path}" | \
        paste-col.py --header "VariantCaller" -v "${caller}" | \
        paste-col.py --header "System" -v "${localhostname}" > \
        "${prefix}.reformat.tsv"
        """
    else if( caller == 'MuTect2' )
        """
        reformat-vcf-table.py -c MuTect2 -s "${sampleID}" -i "${sample_tsv}" | \
        paste-col.py --header "Sample" -v "${sampleID}"  | \
        paste-col.py --header "Run" -v "${params.runID}" | \
        paste-col.py --header "Results" -v "${resultsID}" | \
        paste-col.py --header "Location" -v "${current_dir_path}" | \
        paste-col.py --header "VariantCaller" -v "${caller}" | \
        paste-col.py --header "System" -v "${localhostname}" > \
        "${prefix}.reformat.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process annotate_vcf {
    // annotate the VCF file
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/annotate_vcf", overwrite: true
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(sample_tsv), file(annovar_db_dir) from vcfs_tsvs_reformat.combine(annovar_db_dir)

    output:
    set val(caller), val(sampleID), file(sample_vcf), file(sample_tsv), file(annovar_output_txt), file("${prefix}.avinput.tsv") into vcfs_tsvs_annotations

    script:
    prefix = "${sampleID}.${caller}"
    avinput_file = "${prefix}.avinput"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annovar_output_vcf = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    if( caller == 'HaplotypeCaller' )
        """
        # make sure there are variants present, by checking the .TSV file; should have >1 line
        [ ! "\$( cat "${sample_tsv}" | wc -l )" -gt 1 ] && echo "ERROR: No variants present in ${sample_tsv}, skipping annotation..." && exit 11 || :

        # annovate
        table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --onetranscript \
        --outfile "${prefix}"

        printf "Chr\tStart\tEnd\tRef\tAlt\tAF\tQuality\tAD.ALT\tCHROM\tPOS\tID\tREF\tALT\tQUAL\n" > "${prefix}.avinput.tsv"
        cut -f1-14 ${avinput_file} >>  "${prefix}.avinput.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        # make sure there are variants present, by checking the .TSV file; should have >1 line
        [ ! "\$( cat "${sample_tsv}" | wc -l )" -gt 1 ] && echo "ERROR: No variants present in ${sample_tsv}, skipping annotation..." && exit 11 || :

        table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --onetranscript \
        --outfile "${prefix}"

        printf "Chr\tStart\tEnd\tRef\tAlt\tId\tQuality\tDP\tCHROM\tPOS\tID\tREF\tALT\tQUAL\n" > "${prefix}.avinput.tsv"
        cut -f1-14 ${avinput_file} >>  "${prefix}.avinput.tsv"
        """
    else if( caller == 'MuTect2' )
        """
        # make sure there are variants present, by checking the .TSV file; should have >1 line
        [ ! "\$( cat "${sample_tsv}" | wc -l )" -gt 1 ] && echo "ERROR: No variants present in ${sample_tsv}, skipping annotation..." && exit 11 || :

        table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --onetranscript \
        --outfile "${prefix}"

        # TODO: Need to check this! Need a MuTect2 .vcf with passing variants!
        printf "Chr\tStart\tEnd\tRef\tAlt\tId\tQuality\tDP\tCHROM\tPOS\tID\tREF\tALT\tQUAL\n" > "${prefix}.avinput.tsv"
        cut -f1-14 ${avinput_file} >>  "${prefix}.avinput.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process merge_tables {
    // merge the annotation and vcf tables
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/merge_tables", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(sample_tsv), file(annovar_txt), file(avinput_file) from vcfs_tsvs_annotations

    output:
    set val(caller), val(sampleID), file("${prefix}.annotations.tsv") into merged_tables
    file("${prefix}.annotations.tsv") into merged_tables2

    script:
    prefix = "${sampleID}.${caller}"
    """
    merge-vcf-tables.R "${sample_tsv}" "${annovar_txt}" "${avinput_file}" "${prefix}.annotations.tmp"
    hash-col.py -i "${prefix}.annotations.tmp" -o "${prefix}.annotations.tsv" --header 'Hash' -k Chr Start End Ref Alt CHROM POS REF ALT Sample Run Results VariantCaller
    """

}

process tsv_2_sqlite {
    // convert TSV files into SQLite databases
    // NOTE: case-insensitive columns; 'Ref' and 'REF', one will get dropped...
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/samples/${sampleID}/${caller}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/analysis/tsv_2_sqlite", overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from merged_tables

    output:
    set val(caller), val(sampleID), file("${sqlite}") into samples_sqlite
    file("${dump_csv}")
    file("${dump_sqlite}")

    script:
    prefix = "${sampleID}.${caller}"
    sqlite = "${prefix}.sqlite"
    dump_csv = "${prefix}.sqlite.csv"
    dump_sqlite = "${prefix}.sqlite.txt"
    """
    table2sqlite.py -i "${sample_tsv}" -o "${sqlite}" -t variants --dump-csv "${dump_csv}" --dump-sqlite "${dump_sqlite}"
    """
}

process collect_annotation_tables {
    publishDir "${params.output_dir}/analysis/collect_annotation_tables", mode: 'copy', overwrite: true

    input:
    file('table*') from merged_tables2.collect()

    output:
    file('all_annotations.tsv')

    script:
    """
    concat-tables.py * > all_annotations.tsv
    """
}
