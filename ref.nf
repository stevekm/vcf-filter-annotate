Channel.from("foo").set { input_channel }

process make_ref {
    echo true
    storeDir "${params.ref_dir}"

    input:
    val(foo) from input_channel

    output:
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict")
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")

    script:
    """
    wget https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref.tar.gz
	tar -xvzf ref.tar.gz ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
	rm -f ref.tar.gz
    mv ref/iGenomes .
    ls -l
    pwd
    """
}
