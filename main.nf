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
        SDRF_FOR_COUNT
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

process get_barcodes_filename {
    
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
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile from FASTQ_CDNA_BARCODES_FILES

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

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
        .groupTuple()
        .map{ row-> tuple( row[0], row[1][0]) }
        .set{ TECHREPS }

    // The target set of results will now be the technical replicate group number

    SDRF_FOR_COUNT
        .map{ row-> tuple(row["${params.fields.techrep}"]) }
        .unique()
        .count()
        .set { TARGET_RESULT_COUNT }
    
    // Now add the tech rep group to the run info, group by it, and create a
    // tuple of files keyed by techrep group

    TECHREPS.join( DOWNLOADED_FASTQS )
        .groupTuple(by: 1)
        .map{ row-> tuple( row[1], row[2].flatten(), row[3].flatten()) }
        .set{
            FINAL_FASTQS
        }
}else{
    DOWNLOADED_FASTQS.set{ FINAL_FASTQS }
    
    SDRF_FOR_COUNT
      .map{ row-> tuple(row["${params.fields.run}"]) }
      .unique()
      .count()
      .set { TARGET_RESULT_COUNT }
}

// Remove anything from the cDNA that's not present in the GTF. Otherwise our
// transcript_to_gene mapping will not contain all transcript IDs, and Alevin
// gets upset. This also generates our transcript/ gene mappings

process synchronise_cdna_gtf {

    conda "${baseDir}/envs/cdna_gtf.yml"

    cache 'deep'

    memory { 5.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        file(referenceFasta) from REFERENCE_FASTA.first()
        file(referenceGtf) from REFERENCE_GTF.first()

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

    cache 'deep'
    
    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(referenceFasta) from REFERENCE_FASTA_CLEANED

    output:
        file('salmon_index') into SALMON_INDEX

    """
    salmon index -t ${referenceFasta} -i salmon_index -k ${params.salmon.index.kmerSize}
    """
}

// Run Alevin per row

process alevin {

    conda "${baseDir}/envs/alevin.yml"

    publishDir "$resultsRoot/alevin", mode: 'copy', overwrite: true
    
    memory { 20.GB * task.attempt }
    cpus 12

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(indexDir) from SALMON_INDEX
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz") from FINAL_FASTQS
        file(transcriptToGene) from TRANSCRIPT_TO_GENE

    output:
        set val(runId), file("${runId}") into ALEVIN_RESULTS

    script:

        def alevinType = ''

        if ( params.containsKey(protocol) ){
            alevinType = params.get(protocol).alevinType
        }

    """
    salmon alevin -l ${params.salmon.libType} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        --${alevinType} -i ${indexDir} -p ${task.cpus} -o ${runId} --tgMap ${transcriptToGene}
    """
}

// Convert Alevin output to MTX

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin.yml"
    
    memory { 40.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        set val(runId), file(alevinResults) from ALEVIN_RESULTS

    output:
        file("${runId}_alevin_mtx") into ALEVIN_RESULTS_MTX

    """
    alevinToMtx.py ${runId} ${runId}_alevin_mtx
    """ 
        
} 

// Check the total number of runs we have 

ALEVIN_RESULTS_MTX
    .count()
    .set{ ALEVIN_RESULTS_COUNT } 

process validate_results {
    
    executor 'local'
    
    input:
        val(kallistoResultCount) from ALEVIN_RESULTS_COUNT 
        val(targetCount) from TARGET_RESULT_COUNT

    output:
        stdout DONE

    """
    if [ "$kallistoResultCount" -ne "$targetCount" ]; then
        echo "Alevin results count of $kallistoResultCount does not match expected results number ($targetCount)" 1>&2
        exit 1
    else
        echo "Alevin results count of $kallistoResultCount matches expected results number ($targetCount)"
    fi
    """
}   
