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
            ['LoFreq', 'SC-SERACARE', file("input/SC-SERACARE.LoFreq.vcf.gz")]
            // [ "LoFreq", "NTC-H2O", file("input/vcf_lofreq/NTC-H2O.vcf.gz") ],
            // [ "LoFreq", "HapMap-B17-1267", file("input/vcf_lofreq/HapMap-B17-1267.vcf.gz") ],
            // [ "LoFreq", "SeraCare-1to1", file("input/vcf_lofreq/SeraCare-1to1-Positive.vcf.gz") ],
            // [ "HaplotypeCaller", "NTC-H2O", file("input/vcf_hc/NTC-H2O.vcf.gz") ],
            // [ "HaplotypeCaller", "HapMap-B17-1267", file("input/vcf_hc/HapMap-B17-1267.vcf.gz") ],
            // [ "HaplotypeCaller", "SeraCare-1to1", file("input/vcf_hc/SeraCare-1to1-Positive.vcf.gz") ]
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
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf) from sample_variants_zipped.concat(pairs_variants_zipped)

    output:
    set val(caller), val(sampleID), file("${sampleID}.vcf") into sample_variants

    script:
    """
    gunzip -c "${sample_vcf}" > "${sampleID}.vcf"
    """
}

process normalize_vcf {
    // normalize and split the VCF entries
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta) from sample_variants.combine(ref_fasta)

    output:
    set val(caller), val(sampleID), file("${sampleID}.norm.vcf") into (normalized_variants, normalized_variants2)
    file("${sampleID}.bcftools.multiallelics.stats.txt")
    file("${sampleID}.bcftools.realign.stats.txt")

    script:
    """
    cat ${sample_vcf} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${sampleID}.bcftools.multiallelics.stats.txt" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${sampleID}.bcftools.realign.stats.txt" > \
    "${sampleID}.norm.vcf"
    """
}

process filter_vcf {
    // filter the VCF entries
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true
    echo true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from normalized_variants2.combine(ref_fasta2).combine(ref_fai2).combine(ref_dict2)

    output:
    set val(caller), val(sampleID), file("${sampleID}.filtered.vcf") into (filtered_vcfs, filtered_vcfs2)

    script:
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
            > "${sampleID}.filtered.vcf"
        """
    else if( caller == 'LoFreq' )
        """
        # do not report if frequency is less than 1%
        gatk.sh -T SelectVariants \
                -R "${ref_fasta}" \
                -V "${sample_vcf}" \
                -select "AF > 0.01"  \
                > "${sampleID}.filtered.vcf"
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
                > "${sampleID}.filtered.vcf"
        """
    else
        error "Invalid caller: ${caller}"
}

process vcf_2_tsv {
    // convert VCF file to TSV
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcfs.combine(ref_fasta3).combine(ref_fai3).combine(ref_dict3)

    output:
    set val(caller), val(sampleID), file("${sampleID}.tsv") into vcf_tsvs

    script:
    if( caller == 'HaplotypeCaller' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
        -GF AD -GF DP \
        -o "${sampleID}.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F SB -F INDEL -F CONSVAR -F HRUN \
        -o "${sampleID}.tsv"
        """
    else if( caller == 'MuTect2' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
        -GF AD -GF DP -GF AF \
        -o "${sampleID}.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process reformat_vcf_tsv {
    // reformat and adjust the TSV table for consistency downstream
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from vcf_tsvs

    output:
    set val(caller), val(sampleID), file("${sampleID}.reformat.tsv") into samples_recalc_tsvs

    script:
    if( caller == 'HaplotypeCaller' )
        """
        reformat-vcf-table.py -c GATKHC \
        -s "${sampleID}" \
        -i "${sample_tsv}" \
        -o "${sampleID}.reformat.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        reformat-vcf-table.py -c LoFreq \
        -s "${sampleID}" \
        -i "${sample_tsv}" \
        -o "${sampleID}.reformat.tsv"
        """
    else if( caller == 'MuTect2' )
        """
        reformat-vcf-table.py -c MuTect2 \
        -s "${sampleID}" \
        -i "${sample_tsv}" \
        -o "${sampleID}.reformat.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process update_tsv_keys {
    // add extra columns to the VCF TSV file for downstream
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from samples_recalc_tsvs

    output:
    set val(caller), val(sampleID), file("${sampleID}.updated.tsv") into samples_updated_tsvs

    script:
    """
    # add a column with the sample ID
    paste-col.py -i "${sample_tsv}" \
    -o "${sampleID}.updated.tmp1" \
    --header "Sample" \
    -v "${sampleID}" \
    -d "\t"

    # add a column with the run ID
    paste-col.py -i "${sampleID}.updated.tmp1" \
    -o "${sampleID}.updated.tmp2" \
    --header "Run" \
    -v "${params.runID}" \
    -d "\t"

    # add a column with the results ID
    paste-col.py -i "${sampleID}.updated.tmp2" \
    -o "${sampleID}.updated.tmp3" \
    --header "Results" \
    -v "${resultsID}" \
    -d "\t"

    # add a column with the current dir
    paste-col.py -i "${sampleID}.updated.tmp3" \
    -o "${sampleID}.updated.tmp4" \
    --header "Location" \
    -v "${current_dir_path}" \
    -d "\t"

    # add a column with the variant caller
    paste-col.py -i "${sampleID}.updated.tmp4" \
    -o "${sampleID}.updated.tmp5" \
    --header "VariantCaller" \
    -v "${caller}" \
    -d "\t"

    # add a column with the system hostname
    paste-col.py -i "${sampleID}.updated.tmp5" \
    -o "${sampleID}.updated.tmp6" \
    --header "System" \
    -v "${localhostname}" \
    -d "\t"

    # add md5sum has for desired columns
    hash-col.py -i "${sampleID}.updated.tmp6" \
    -o "${sampleID}.updated.tsv" \
    --header 'Hash' \
    -k CHROM POS REF ALT Sample Run Results VariantCaller
    """
}

process annotate_vcf {
    // annotate the VCF file
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(annovar_db_dir) from filtered_vcfs2.combine(annovar_db_dir)

    output:
    set val(caller), val(sampleID), file(annovar_output_txt), file(avinput_file) into samples_annotations
    file(annovar_output_vcf)
    set val(caller), val(sampleID), file(avinput_file) into samples_avinput

    script:
    avinput_file = "${sampleID}.avinput"
    annovar_output_txt = "${sampleID}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annovar_output_vcf = "${sampleID}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    """
    # convert to ANNOVAR format
    # convert2annovar.pl --format vcf4old --includeinfo --comment "${sample_vcf}" --outfile "${avinput_file}"
    # check number of lines between the files
    # [ ! "\$( cat "${avinput_file}" | wc -l )" -eq "\$(grep -v '^#' "${sample_vcf}" | wc -l)" ] && echo "ERROR: number of entries does not match between files ${sample_vcf} and ${avinput_file}" && exit 1 || :

    # annovate
    table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
    --buildver "${params.ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${params.ANNOVAR_PROTOCOL}" \
    --operation "${params.ANNOVAR_OPERATION}" \
    --nastring . \
    --vcfinput \
    --onetranscript \
    --outfile "${sampleID}"
    """
    // --otherinfo \
}

process reformat_avinput {
    // add column headers for the .avinput file and subset it for just the desired columns
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(avinput_file) from samples_avinput

    when:
    caller == 'HaplotypeCaller'

    output:
    set val(caller), val(sampleID), file("${sampleID}.avinput.reformat.tsv")

    script:
    if( caller == 'HaplotypeCaller' )
        """
        # keep only the first columns and add headers
        printf "Chr\tStart\tEnd\tRef\tAlt\tAF\tQUAL\tAD.ALT\tCHROM\tPOS\tID\tREF\tALT\n" > "${sampleID}.avinput.reformat.tsv"
        cut -f1-13 ${avinput_file} >>  "${sampleID}.avinput.reformat.tsv"
        """
}

process tsv_2_sqlite {
    // convert TSV files into SQLite databases
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from samples_updated_tsvs

    output:
    set val(caller), val(sampleID), file("${sampleID}.sqlite") into samples_sqlite

    script:
    """
    table2sqlite.py -i "${sample_tsv}" -o "${sampleID}.sqlite" -t variants 
    """
}

// samples_annotations.join(samples_recalc_tsvs, by: [0,1]).tap { samples_annotations_tables }
// process merge_annotation_tables {
//     tag "${caller}-${sampleID}"
//     publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true
//
//     input:
//     set val(caller), val(sampleID), file(annovar_output_txt), file(avinput_file), file(recalc_tsv) from samples_annotations_tables
//
//     output:
//     file(output_file)
//
//     when:
//     caller == 'HaplotypeCaller'
//
//     script:
//     output_file = "${sampleID}.annotation.tsv"
//     if( caller == 'HaplotypeCaller' )
//         """
//         merge-ANNOVAR-tables-GATKHC.R "${sampleID}" "${recalc_tsv}" "${annovar_output_txt}" "${avinput_file}" "${output_file}"
//         """
//     else
//         error "Invalid caller: ${caller}"
// }
