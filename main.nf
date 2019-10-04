#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
referenceGtf = params.referenceGtf
protocol = params.protocol

manualDownloadFolder =''
if ( params.containsKey('manualDownloadFolder')){
    manualDownloadFolder = params.manualDownloadFolder
}

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

REFERENCE_FASTA = Channel.fromPath( referenceFasta, checkIfExists: true )
REFERENCE_GTF = Channel.fromPath( referenceGtf, checkIfExists: true )

// Read URIs from SDRF, generate target file names, and barcode locations

SDRF_FOR_FASTQS
    .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.cdna_uri}"], row["${params.fields.cell_barcode_uri}"], file(row["${params.fields.cdna_uri}"]).getName(), file(row["${params.fields.cell_barcode_uri}"]).getName(), row["${params.fields.cell_barcode_size}"], row["${params.fields.umi_barcode_size}"], row["${params.fields.end}"], row["${params.fields.cell_count}"]) }
    .set { FASTQ_RUNS }

// Call the download script to retrieve run fastqs

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 1.hour * task.attempt }

    errorStrategy { task.attempt<=10 ? 'retry' : 'finish' } 
    
    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile, val(barcodeLength), val(umiLength), val(end), val(cellCount) from FASTQ_RUNS

    output:
        set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}"), val(barcodeLength), val(umiLength), val(end), val(cellCount) into DOWNLOADED_FASTQS

    """
        if [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
           ln -s $manualDownloadFolder/${cdnaFastqFile} ${cdnaFastqFile}
           ln -s $manualDownloadFolder/${barcodesFastqFile} ${barcodesFastqFile}
        elif [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ ! -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is available locally, but barcodes file $barcodesFastqFile is not 1>&2
            exit 2    
        elif [ -n "$manualDownloadFolder" ] && [ ! -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is not available locally, but barcodes file $barcodesFastqFile is 1>&2
            exit 3    
        else
            confPart=''
            if [ -e "$NXF_TEMP/atlas-fastq-provider/download_config.sh" ]; then
                confPart=" -c $NXF_TEMP/atlas-fastq-provider/download_config.sh"
            fi 
            fetchFastq.sh -f ${cdnaFastqURI} -t ${cdnaFastqFile} -m ${params.downloadMethod} \$confPart
            fetchFastq.sh -f ${barcodesFastqURI} -t ${barcodesFastqFile} -m ${params.downloadMethod} \$confPart
        fi
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
        .map{ row-> tuple( row[1], row[2].flatten(), row[3].flatten(), row[4][0], row[5][0], row[6][0], row[7][0]) }
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

REFERENCE_FASTA_CLEANED
    .into{
        REFERENCE_FASTA_CLEANED_FOR_SALMON
        REFERENCE_FASTA_CLEANED_FOR_KALLISTO
    }

// Generate an index from the transcriptome

process salmon_index {

    conda "${baseDir}/envs/alevin.yml"

    cache 'deep'
    
    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(referenceFasta) from REFERENCE_FASTA_CLEANED_FOR_SALMON

    output:
        file('salmon_index') into SALMON_INDEX

    """
    salmon index -t ${referenceFasta} -i salmon_index -k ${params.salmon.index.kmerSize}
    """
}

process kallisto_index {
    
    conda "${baseDir}/envs/kallisto.yml"

    cache 'deep'
    
    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(referenceFasta) from REFERENCE_FASTA_CLEANED_FOR_KALLISTO

    output:
        file('kallisto_index.idx') into KALLISTO_INDEX

    """
    kallisto index -i kallisto_index.idx -k ${params.salmon.index.kmerSize} ${referenceFasta}
    """
}

FINAL_FASTQS
    .into{
        FASTQS_FOR_ALEVIN
        FASTQS_FOR_KALLISTO
    }


// Run Alevin per row

process alevin {

    conda "${baseDir}/envs/alevin.yml"
    
    cache 'deep'

    memory { 20.GB * task.attempt }
    cpus 12

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(indexDir) from SALMON_INDEX
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FASTQS_FOR_ALEVIN
        file(transcriptToGene) from TRANSCRIPT_TO_GENE

    output:
        set val('alevin'), val(runId), file("${runId}"),  file("${runId}_pre/alevin/raw_cb_frequency.txt"), file("${runId}/alevin/quants_mat.gz"), file("${runId}/alevin/quants_mat_cols.txt"), file("${runId}/alevin/quants_mat_rows.txt") into ALEVIN_RESULTS

    script:

        def barcodeConfig = ''

        if ( params.containsKey(protocol) ){

            canonicalProtocol = params.get(protocol)
            alevinType = canonicalProtocol.alevinType

            // Non-standard barcode config is supplied as a custom method

            if ( "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
                barcodeConfig = "--barcodeLength ${barcodeLength} --umiLength ${umiLength} --end ${end}" 

            }else{
                barcodeConfig = "--$alevinType"
            }
            
        }

    """
    if [ -z "$barcodeConfig" ]; then
        echo Input of $protocol results is misconfigured 1>&2
        exit 1
    fi

    # Do a pre-run to derive a starting whitelist, see https://github.com/COMBINE-lab/salmon/issues/362

    salmon alevin -l ${params.salmon.libType} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        ${barcodeConfig} -i ${indexDir} -p ${task.cpus} -o ${runId}_pre --tgMap ${transcriptToGene} --dumpFeatures --noQuant
   
    # Derive a relaxed whitelist, removing only the most obvious junk 

    if [ \$? -eq 0 ]; then 
        awk '{ if (\$2 > ${params.minCbFreq }) { print \$1} }' ${runId}_pre/alevin/raw_cb_frequency.txt > pre_whitelist.txt
    fi

    # Supply the whitelist to the main Alevin run

    salmon alevin -l ${params.salmon.libType} -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        ${barcodeConfig} -i ${indexDir} -p ${task.cpus} -o ${runId}_tmp --tgMap ${transcriptToGene} --whitelist pre_whitelist.txt \
        --forceCells \$(cat pre_whitelist.txt | wc -l | tr -d '\\n') --dumpMtx
 
    mv ${runId}_tmp ${runId}
    """
}

process kallisto_bus {
    
    conda "${baseDir}/envs/kallisto.yml"
    
    cache 'deep'

    memory { 20.GB * task.attempt }
    cpus 12

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        file(indexFile) from KALLISTO_INDEX
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FASTQS_FOR_KALLISTO
        file(transcriptToGene) from TRANSCRIPT_TO_GENE

    output:
        set val(runId), file("${runId}") into KALLISTO_RESULTS

    script:
        
        def barcodeConfig = ''

        if ( params.containsKey(protocol) ){

            canonicalProtocol = params.get(protocol)
            kallistoType = canonicalProtocol.kallistoType
            umiEnd = barcodeLength + umiLength

            // Non-standard barcode config is supplied as a custom method

            if ( "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
                barcodeConfig = "-x 0,0,$barcodeLength:0,$barcodeLength,$umiLength:1,0,0" 

            }else{
                barcodeConfig = "-x $kallistoType"
            }
            
        }

        """
        mkdir -p ${runId}_tmp

        kallisto bus -i ${indexFile} -o ${runId}_tmp $barcodeConfig -t ${task.cpus} \$(paste -d ' ' <(ls barcodes*.fastq.gz) <(ls cdna*.fastq.gz) | tr '\\n' ' ')
        mkdir -p ${runId}_tmp/counts_mtx
        bustools sort -T \$(pwd)/tmp/ -t ${task.cpus} -p ${runId}_tmp/output.bus | \
            bustools count -o ${runId}_tmp/counts_mtx/gene -g $transcriptToGene -e ${runId}_tmp/matrix.ec -t ${runId}_tmp/transcripts.txt --genecounts -
        gzip ${runId}_tmp/counts_mtx/gene.mtx        

        mv ${runId}_tmp ${runId}

        """ 
}

// Just a utility process to make a raw cb frequency list for comparision with Alevin's

process add_raw_counts_to_kallisto {
    
    conda "${baseDir}/envs/kallisto.yml"

    input:
        set val(runId), file(kallistoDir) from KALLISTO_RESULTS

    output:
        set val('kallisto'), val(runId), file(kallistoDir), file('raw_cb_frequency.txt'), file("${runId}/counts_mtx/gene.mtx.gz"), file("${runId}/counts_mtx/gene.genes.txt"), file("${runId}/counts_mtx/gene.barcodes.txt") into KALLISTO_RESULTS_WITH_FREQS
    

    """
    bustools text ${runId}/output.bus -p |  awk '{print \$1}' | sort | uniq -c | sort -nr | awk '{print \$2"\t"\$1}' > raw_cb_frequency.txt
    """

}

ALEVIN_RESULTS
    .concat( KALLISTO_RESULTS_WITH_FREQS)
    .into{
        RESULTS_FOR_QC
        RESULTS_FOR_PROCESSING
        RESULTS_FOR_OUTPUT
    }


// Convert Alevin output to MTX. There will be one of these for every run, or
// technical replicate group of runs

process convert_to_mtx {

    conda "${baseDir}/envs/parse_alevin.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(method), val(runId), file(resultDir), file(rawBarcodeFreq), file(mtxFile), file(genesFile), file(barcodesFile) from RESULTS_FOR_PROCESSING

    output:
        set val(method), val(runId), file("counts_mtx") into RAW_MTX

    """
    mtxTo10x.py --cell_prefix ${runId}- $mtxFile $genesFile $barcodesFile  counts_mtx
    """ 
}

RAW_MTX
    .into{
        RAW_MTX_FOR_QC
        RAW_MTX_FOR_EMPTYDROPS
        RAW_MTX_FOR_OUTPUT
    }

// Make a diagnostic plot

RESULTS_FOR_QC
    .join(RAW_MTX_FOR_QC, by: [0,1])
    .set{
        QC_INPUTS
    }

process droplet_qc_plot{
    
    conda "${baseDir}/envs/alevin.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(method), val(runId), file(alevinResult), file(rawBarcodeFreq), file(mtx) from QC_INPUTS

    output:
        set val(method), val(runId), file("${method}-${runId}.png") into QC_PLOTS

    """
    dropletBarcodePlot.R $rawBarcodeFreq $mtx $runId ${method}-${runId}.png
    """ 
}

// Remove empty droplets from Alevin results

process remove_empty_drops {
    
    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'ignore' }
    maxRetries 20
   
    input:
        set val(method), val(runId), file(countsMtx) from RAW_MTX_FOR_EMPTYDROPS

    output:
        set val(method), val(runId), file('nonempty.rds') into NONEMPTY_RDS

    """
        dropletutils-read-10x-counts.R -s counts_mtx -c TRUE -o matrix.rds
        dropletutils-empty-drops.R -i matrix.rds --lower ${params.emptyDrops.lower} --niters ${params.emptyDrops.nIters} --filter-empty ${params.emptyDrops.filterEmpty} \
            --filter-fdr ${params.emptyDrops.filterFdr} --ignore ${params.minCbFreq} -o nonempty.rds -t nonempty.txt
    """
}

// Convert R matrix object with filtered cells back to .mtx

process rds_to_mtx{

    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
   
    input:
        set val(method), val(runId), file(rds) from NONEMPTY_RDS

    output:
        set val(method), val(runId), file("counts_mtx_nonempty") into NONEMPTY_MTX

    """ 
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))

        counts_sce <- readRDS('$rds')
        write10xCounts(assays(counts_sce)[[1]], path = 'counts_mtx_nonempty', barcodes = colData(counts_sce)\$Barcode, gene.id = rownames(counts_sce))
    """
}

// Compile raw results with raw and emptyDrops-filtered MTX

RESULTS_FOR_OUTPUT
    .join(RAW_MTX_FOR_OUTPUT, by: [0,1])
    .join(NONEMPTY_MTX, by: [0,1])
    .join(QC_PLOTS, by: [0,1])
    .set{ COMPILED_RESULTS }

process compile_results{

    publishDir "$resultsRoot/alevin", mode: 'copy', overwrite: true
    
    input:
        set val(method), val(runId), file('raw_alevin'), file(rawBarcodeFreq), file(countsMtx), file(countsMtxNonempty), file(qcPlot) from COMPILED_RESULTS

    output:
        set val(method), val(runId), file("$runId") into RESULTS_FOR_COUNTING

    """
        mkdir -p raw_$method/$method/mtx
        cp -P $countsMtx $countsMtxNonempty raw_$method/$method/mtx 
        mkdir -p raw_$method/$method/qc
        cp -P $qcPlot raw_$method/$method/qc
        cp -P raw_$method $runId
    """
}

// Check the total number of runs we have 

RESULTS_FOR_COUNTING
    .count()
    .set{ RESULTS_COUNT } 

process validate_results {
    
    executor 'local'
    
    input:
        val(resultCount) from RESULTS_COUNT 
        val(targetCount) from TARGET_RESULT_COUNT

    output:
        stdout DONE

    """
    expectedCount=\$(( 2*$targetCount )) 

    if [ "\$expectedCount" -ne "$targetCount" ]; then
        echo "Alevin and Kallisto results count of $resultCount does not match expected results number (\$expectedCount)" 1>&2
        exit 1
    else
        echo "Alevin and Kallisto results count of $resultCount matches expected results number (\$expectedCount)"
    fi
    """
}   

