#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
referenceGtf = params.referenceGtf
protocol = params.protocol

// Read ENA_RUN column from an SDRF

Channel
    .fromPath(sdrfFile, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .filter{ row -> (! row.containsKey(params.fields.quality)) || ( row["${params.fields.quality}"].toLowerCase() != 'not ok') }
    .into {
        SDRF_FOR_FASTQS
        SDRF_FOR_STRAND
        SDRF_FOR_TECHREP
    }

SDRF_FOR_FASTQS
    .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.cdna_uri}"], row["${params.fields.cell_barcode_uri}"]) }
    .set { FASTQ_RUNS }

REFERENCE_FASTA = Channel.fromPath( referenceFasta, checkIfExists: true )
REFERENCE_GTF = Channel.fromPath( referenceGtf, checkIfExists: true )

// Get the file names from the URLs

process get_cdna_filename {
    
    executor 'local'

    input:
        set runId, cdnaFastqURI, barcodesFastqURI from FASTQ_RUNS
    
    output:
        set val(runId), val(cdnaFastqURI), val(barcodesFastqURI), stdout into FASTQ_CDNA_FILES

    """
        basename $cdnaFastqURI | tr -d \'\\n\'
    """
}

process get_cdna_filename {
    
    executor 'local'

    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile from FASTQ_CDNA_FILES
    
    output:
        set runId, val(cdnaFastqURI), val(barcodesFastqURI), val(cdnaFastqFile), stdout into FASTQ_CDNA_BARCODES_FILES

    """
        basename $barcodesFastqURI | tr -d \'\\n\'
    """
}

// Call the download script to retrieve run fastqs

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 1.hour * task.attempt }

    errorStrategy { task.attempt<=10 ? 'retry' : 'finish' } 
    
    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile into FASTQ_CDNA_BARCODES_FILES

    output:
        set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}") into DOWNLOADED_FASTQS

    """
        confPart=''
        if [ -e "$NXF_TEMP/atlas-fastq-provider/download_config.sh" ]; then
            confPart=" -c $NXF_TEMP/atlas-fastq-provider/download_config.sh"
        fi 
        fetchFastq.sh -f ${cdnaFastqURI} -t ${cdnaFastqFile} -m ${params.downloadMethod} \$confPart
        fetchFastq.sh -f ${barcodesFastqURI} -t ${barcodesFastqFile} -m ${params.downloadMethod} \$confPart
    """
}

// Remove anything from the cDNA that's not present in the GTF. Otherwise our
// transcript_to_gene mapping will not contain all transcript IDs, and Alevin
// gets upset. This also generates our transcript/ gene mappings

process synchronise_cdna_gtf {

    conda "${baseDir}/envs/cdna_gtf.yml"

    input:
        file(referenceFasta) from REFERENCE_FASTA
        file(referenceGtf) from REFERENCE_GTF

    output:
        file('transcript_to_gene.txt') into TRANSCRIPT_TO_GENE
        file('cleanedCdna.fa.gz') into REFERENCE_FASTA_CLEANED

    """
    setupCdnaTranscriptToGene.R ${referenceFasta} ${referenceGtf} transcript_to_gene.txt cleanedCdna.fa.gz 
    """

}

// Generate an index from the transcriptome

process salmon_index {

    conda "${baseDir}/envs/alevin.yml"

    input:
        file(referenceFasta) from REFERENCE_FASTA_CLEANED

    output:
        file('salmon_index')

    """
    salmon index -t ${referenceFasta} -i salmon_index -k ${params.index.kmer_size}
    """
}

// Run Alevin per row

process alevin {

    conda "${baseDir}/envs/alevin.yml"

    input:
        file(indexDir) from SALMON_INDEX
        set val(runId), file(cdnaFastqFile), file(barcodesFastqFile) from DOWNLOADED_FASTQS
        file(transcriptToGene) from TRANSCRIPT_TO_GENE

    output:
        set val(runId), file("${runId}_alevin") into ALEVIN_RESULTS

    script:

        def alevinType = ''

        if ( params.containsKey(protocol) ){
            alevinType = params.get(protocol).alevinType
        }

    """
    salmon alevin -l ISR -1 ${barcodesFastqFile} -2 ${cdnaFastqFile} --${alevinType} -i ${indexDir} -p ${task.cpus} \
        -o ${runId}_alevin --tgMap ${transcriptToGene}
    """
}

# Convert Alevin output to MTX

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin.yml"

    input:
        set val(runId), file(alevinResults) from ALEVIN_RESULTS

    output:
        file("${run_id}_alevin_mtx") into ALEVIN_RESULTS_MTX

    """
    alevinToMtx.py ${runId}_alevin ${runId}_alevin_mtx
    """ 
        
} 
