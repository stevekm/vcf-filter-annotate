Channel.from( [ ['HaplotypeCaller', 'NC-HAPMAP', file("input/NC-HAPMAP.original.vcf")] ] )
        .set { sample_variants }
Channel.fromPath("ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
        .into { ref_fasta; ref_fasta2; ref_fasta3 }

Channel.fromPath("ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")
        .into { ref_fai; ref_fai2; ref_fai3 }
Channel.fromPath("ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict")
        .into { ref_dict; ref_dict2; ref_dict3 }
Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }

// # 1. left normalize indels & split multiallelic entries (.vcf -> norm.vcf)
// # 2. filtere vcf (norm.vcf -> .norm.filtered.vcf)
// # 2. convert to tsv (.norm.filtered.vcf -> .tsv)
// # 3. recalculate allele frequency (.tsv -> recalc.tsv)
// # 4. annotate (norm.vcf -> .avinput, .hg19_multianno.txt)
// # 5. merge annotations with .tsv (.hg19_multianno.txt, recalc.tsv -> ... )

process normalize_vcf {
    tag "${sampleID}"
    publishDir "${params.output_dir}/normalize_vcf", mode: 'copy', overwrite: true

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
    tag "${sampleID}"
    echo true

    input:
    set val(caller), val(sampleID), file(sample_vcf) from normalized_variants

    script:
    """
    grep '6676635' "${sample_vcf}"
    """
}


process filter_vcf {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/filter_vcf", mode: 'copy', overwrite: true
    echo true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from normalized_variants2.combine(ref_fasta2).combine(ref_fai2).combine(ref_dict2)

    output:
    set val(caller), val(sampleID), file("${sampleID}.norm.filtered.vcf") into (filtered_vcfs, filtered_vcfs2)

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
            > "${sampleID}.norm.filtered.vcf"
        """
    else
        error "Invalid caller: ${caller}"
}

process vcf_2_tsv {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/vcf_2_tsv", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcfs.combine(ref_fasta3).combine(ref_fai3).combine(ref_dict3)

    output:
    set val(caller), val(sampleID), file("${sampleID}.norm.filtered.tsv") into vcf_tsvs

    script:
    if( caller == 'HaplotypeCaller' )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${sample_vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
        -GF AD -GF DP \
        -o "${sampleID}.norm.filtered.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process recalc_tsv {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/recalc_tsv", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_tsv) from vcf_tsvs

    output:
    set val(caller), val(sampleID), file("${sampleID}.norm.filtered.recalc.tsv")

    script:
    if( caller == 'HaplotypeCaller' )
        """
        recalc-vcf-AF.py -c GATKHC -s "${sampleID}" -i "${sample_tsv}" -o "${sampleID}.norm.filtered.recalc.tsv"
        """
    else
        error "Invalid caller: ${caller}"
}

process annotate_vcf {
    tag "${caller}-${sampleID}"
    publishDir "${params.output_dir}/annotate_vcf", mode: 'copy', overwrite: true

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(annovar_db_dir) from filtered_vcfs2.combine(annovar_db_dir)

    output:
    file(annovar_output_file)

    script:
    avinput_file = "${sampleID}.avinput"
    annovar_output_file = "${sampleID}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    """
    # convert to ANNOVAR format
    convert2annovar.pl --format vcf4old --includeinfo "${sample_vcf}" --outfile "${avinput_file}"


    # check number of lines between the files
    [ ! "\$( cat "${avinput_file}" | wc -l )" -eq "\$(grep -v '^#' "${sample_vcf}" | wc -l)" ] && echo "ERROR: number of entries does not match between files ${sample_vcf} and ${avinput_file}" && exit 1 || :

    # annovate
    table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
    --buildver "${params.ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${params.ANNOVAR_PROTOCOL}" \
    --operation "${params.ANNOVAR_OPERATION}" \
    --nastring . \
    --outfile "${sampleID}"
    """
}
