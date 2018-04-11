// ~~~~~  PIPELINE FOR UNPAIRED FILES ~~~~~ //
Channel.from( [
            ['HaplotypeCaller', 'NC-HAPMAP', file("input/NC-HAPMAP.GATKHC.vcf.gz")],
            ['LoFreq', 'NC-HAPMAP', file("input/NC-HAPMAP.LoFreq.vcf.gz")],
            ['HaplotypeCaller', 'SC-SERACARE', file("input/SC-SERACARE.GATKHC.vcf.gz")],
            ['LoFreq', 'SC-SERACARE', file("input/SC-SERACARE.LoFreq.vcf.gz")],
            [ "LoFreq", "NTC-H2O", file("input/vcf_lofreq/NTC-H2O.vcf.gz") ],
            [ "LoFreq", "HapMap-B17-1267", file("input/vcf_lofreq/HapMap-B17-1267.vcf.gz") ],
            [ "LoFreq", "SeraCare-1to1", file("input/vcf_lofreq/SeraCare-1to1-Positive.vcf.gz") ],
            [ "HaplotypeCaller", "NTC-H2O", file("input/vcf_hc/NTC-H2O.vcf.gz") ],
            [ "HaplotypeCaller", "HapMap-B17-1267", file("input/vcf_hc/HapMap-B17-1267.vcf.gz") ],
            [ "HaplotypeCaller", "SeraCare-1to1", file("input/vcf_hc/SeraCare-1to1-Positive.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr8", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr8.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr2", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr2.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr5", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr5.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr12", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr12.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr3", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr3.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr10", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr10.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr20", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr20.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr14", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr14.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chrX", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrX.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr13", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr13.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr19", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr19.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr22", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr22.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr15", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr15.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr9", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr9.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr4", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr4.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr7", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr7.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr16", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr16.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr17", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr17.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr18", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr18.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr11", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr11.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr1", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr1.vcf.gz") ],
            [ "MuTect2","SeraCare-1to1.HapMap-B17-1267.chrY", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrY.vcf.gz") ],
            [ "MuTect2", "SeraCare-1to1.HapMap-B17-1267.chr6", file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr6.vcf.gz") ]
            ] )
        .set { sample_variants_zipped }

Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
        .into { ref_fasta; ref_fasta2; ref_fasta3 }

Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")
        .into { ref_fai; ref_fai2; ref_fai3 }
Channel.fromPath("${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict")
        .into { ref_dict; ref_dict2; ref_dict3 }
Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }

// # 1. left normalize indels & split multiallelic entries (.vcf -> norm.vcf)
// # 2. filtere vcf (norm.vcf -> .norm.filtered.vcf)
// # 2. convert to tsv (.norm.filtered.vcf -> .tsv)
// # 3. recalculate allele frequency (.tsv -> recalc.tsv)
// # 4. annotate (norm.vcf -> .avinput, .hg19_multianno.txt)
// # 5. merge annotations with .tsv (.hg19_multianno.txt, recalc.tsv -> ... )

process unzip_samples {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf) from sample_variants_zipped

    output:
    set val(caller), val(sampleID), file("${sampleID}.vcf") into sample_variants

    script:
    """
    gunzip -c "${sample_vcf}" > "${sampleID}.vcf"
    """
}

process normalize_vcf {
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

process check_normalization {
    tag "${caller}-${sampleID}"
    echo true

    input:
    set val(caller), val(sampleID), file(sample_vcf) from normalized_variants

    script:
    """
    echo "[check_normalization] ['${caller}-${sampleID}'] \$(grep '6676635' '${sample_vcf}')"
    """
}


process filter_vcf {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true
    echo true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from normalized_variants2.combine(ref_fasta2).combine(ref_fai2).combine(ref_dict2)

    output:
    set val(caller), val(sampleID), file("${sampleID}.filtered.vcf") into (filtered_vcfs, filtered_vcfs2)

    when:
    caller != 'MuTect2'

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

            # alternate allele freq (allele depth / depth) greater than 0.5
            # -select "vc.getGenotype('${sampleID}').getAD().1 / vc.getGenotype('${sampleID}').getDP() > 0.50" \
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
    else
        error "Invalid caller: ${caller}"
}

process vcf_2_tsv {
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
    else
        error "Invalid caller: ${caller}"
}

process recalc_tsv {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from vcf_tsvs

    output:
    set val(caller), val(sampleID), file("${sampleID}.recalc.tsv") into samples_recalc_tsvs

    script:
    if( caller == 'HaplotypeCaller' )
        """
        recalc-vcf-AF.py -c GATKHC -s "${sampleID}" -i "${sample_tsv}" -o "${sampleID}.recalc.tsv"
        """
    else if( caller == 'LoFreq' )
        """
        recalc-vcf-AF.py -c LoFreq -s "${sampleID}" -i "${sample_tsv}" -o "${sampleID}.recalc.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process annotate_vcf {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(annovar_db_dir) from filtered_vcfs2.combine(annovar_db_dir)

    output:
    set val(caller), val(sampleID), file(annovar_output_txt), file(avinput_file) into samples_annotations
    file(annovar_output_vcf)

    when:
    caller == 'HaplotypeCaller'

    script:
    avinput_file = "${sampleID}.avinput"
    annovar_output_txt = "${sampleID}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annovar_output_vcf = "${sampleID}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    if( caller == 'HaplotypeCaller' )
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
        --otherinfo \
        --onetranscript \
        --outfile "${sampleID}"
        """
    else
        error "Invalid caller: ${caller}"
}

samples_annotations.join(samples_recalc_tsvs, by: [0,1]).tap { samples_annotations_tables }

process merge_annotation_tables {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(annovar_output_txt), file(avinput_file), file(recalc_tsv) from samples_annotations_tables

    output:
    file(output_file)

    script:
    output_file = "${sampleID}.annotation.tsv"
    if( caller == 'HaplotypeCaller' )
        """
        # merge-ANNOVAR-tables-GATKHC.R sampleID recalc_tsv_file annovar_file avinput_file output_file
        merge-ANNOVAR-tables-GATKHC.R "${sampleID}" "${recalc_tsv}" "${annovar_output_txt}" "${avinput_file}" "${output_file}"
        """
    else
        error "Invalid caller: ${caller}"
}




// ~~~~~  PIPELINE FOR PAIRED FILES ~~~~~ //
// Channel.from([
//     // caller, tumorID, normalID, pairID, file
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr8.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr2.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr5.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr12.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr3.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr10.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr20.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr14.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrX.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr13.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr19.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr22.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr21.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr15.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr9.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr4.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr7.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr16.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr17.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr18.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr11.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr1.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chrY.vcf.gz") ],
// [ "MuTect2", "SeraCare-1to1", "HapMap-B17-1267", "SeraCare-1to1.HapMap-B17-1267" file("input/vcf_mutect2/SeraCare-1to1-Positive_HapMap-B17-1267.chr6.vcf.gz") ]
//     ])
//     .set { pairs_variants_zipped }
//
//
// process unzip_pairs {
//     tag "${caller}-${sampleID}"
//     publishDir "${params.output_dir}/${sampleID}/${caller}", mode: 'copy', overwrite: true
//
//     input:
//     set val(caller), val(sampleID), file(sample_vcf) from sample_variants_zipped
//
//     output:
//     set val(caller), val(sampleID), file("${sampleID}.vcf") into sample_variants
// }
