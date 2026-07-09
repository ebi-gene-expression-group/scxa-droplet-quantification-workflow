#!/usr/bin/env nextflow

WorkflowParamValidator.validate(params)

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
transcriptToGene = params.transcriptToGene
transcriptomeIndex = params.transcriptomeIndex
protocol = params.protocol
experimentType = params.experimentType

manualDownloadFolder =''
if ( params.containsKey('manualDownloadFolder')){
    manualDownloadFolder = params.manualDownloadFolder
}

fastqProviderConfig = ''
if ( params.containsKey('fastqProviderConfig')){
    fastqProviderConfig = params.fastqProviderConfig
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

TRANSCRIPT_TO_GENE = Channel.fromPath( transcriptToGene, checkIfExists: true ).first()

// Read URIs from SDRF, generate target file names, and barcode locations

SDRF_FOR_FASTQS
    .map{ row-> 
      controlled_access = 'no'
      if (  params.fields.containsKey('controlled_access')){
        controlled_access = WorkflowParamValidator.safeControlledAccess(row["${params.fields.controlled_access}"])
      }
      def cdna_uri = WorkflowParamValidator.safeUri(row["${params.fields.cdna_uri}"], params.fields.cdna_uri)
      def cell_barcode_uri = WorkflowParamValidator.safeUri(row["${params.fields.cell_barcode_uri}"], params.fields.cell_barcode_uri)
      tuple(
        WorkflowParamValidator.safeToken(row["${params.fields.run}"], params.fields.run),
        cdna_uri,
        cell_barcode_uri,
        WorkflowParamValidator.safeToken(file(cdna_uri).getName(), "${params.fields.cdna_uri} basename"),
        WorkflowParamValidator.safeToken(file(cell_barcode_uri).getName(), "${params.fields.cell_barcode_uri} basename"),
        WorkflowParamValidator.safeInteger(row["${params.fields.cell_barcode_size}"], params.fields.cell_barcode_size),
        WorkflowParamValidator.safeInteger(row["${params.fields.umi_barcode_size}"], params.fields.umi_barcode_size),
        WorkflowParamValidator.safeInteger(row["${params.fields.end}"], params.fields.end),
        WorkflowParamValidator.safeInteger(row["${params.fields.cell_count}"], params.fields.cell_count),
        controlled_access
      )
    }    
    .set { FASTQ_RUNS }

// Call the download script to retrieve run fastqs

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 10.hour * task.attempt }
    memory { 20.GB * task.attempt }

    errorStrategy { task.attempt<=10 & task.exitStatus != 4 ? 'retry' : 'finish' } 
    
    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile, val(barcodeLength), val(umiLength), val(end), val(cellCount), val(controlledAccess) from FASTQ_RUNS

    output:
        set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}"), val(barcodeLength), val(umiLength), val(end), val(cellCount) into DOWNLOADED_FASTQS

    """
        DOWNLOAD_METHOD=${WorkflowParamValidator.shellQuote(params.downloadMethod)}
        MANUAL_DOWNLOAD_FOLDER=${WorkflowParamValidator.shellQuote(manualDownloadFolder)}
        FASTQ_PROVIDER_CONFIG=${WorkflowParamValidator.shellQuote(fastqProviderConfig)}
        CDNA_FASTQ_FILE=${WorkflowParamValidator.shellQuote(cdnaFastqFile)}
        BARCODES_FASTQ_FILE=${WorkflowParamValidator.shellQuote(barcodesFastqFile)}
        CONTROLLED_ACCESS=${WorkflowParamValidator.shellQuote(controlledAccess)}

        if ! [ -z "$ATLAS_TMPDIR" ]; then 
            TMPDIR=$ATLAS_TMPDIR; 
        else 
            echo "NOTE: ATLAS_TMPDIR not defined"
        fi
        if [ -n "\$MANUAL_DOWNLOAD_FOLDER" ] && [ -e "\$MANUAL_DOWNLOAD_FOLDER/\$CDNA_FASTQ_FILE" ] && [ -e "\$MANUAL_DOWNLOAD_FOLDER/\$BARCODES_FASTQ_FILE" ]; then
           ln -s "\$MANUAL_DOWNLOAD_FOLDER/\$CDNA_FASTQ_FILE" "\$CDNA_FASTQ_FILE"
           ln -s "\$MANUAL_DOWNLOAD_FOLDER/\$BARCODES_FASTQ_FILE" "\$BARCODES_FASTQ_FILE"
        elif [ -n "\$MANUAL_DOWNLOAD_FOLDER" ] && [ -e "\$MANUAL_DOWNLOAD_FOLDER/\$CDNA_FASTQ_FILE" ] && [ ! -e "\$MANUAL_DOWNLOAD_FOLDER/\$BARCODES_FASTQ_FILE" ]; then
            echo "cDNA file \$CDNA_FASTQ_FILE is available locally, but barcodes file \$BARCODES_FASTQ_FILE is not" 1>&2
            exit 2    
        elif [ -n "\$MANUAL_DOWNLOAD_FOLDER" ] && [ ! -e "\$MANUAL_DOWNLOAD_FOLDER/\$CDNA_FASTQ_FILE" ] && [ -e "\$MANUAL_DOWNLOAD_FOLDER/\$BARCODES_FASTQ_FILE" ]; then
            echo "cDNA file \$CDNA_FASTQ_FILE is not available locally, but barcodes file \$BARCODES_FASTQ_FILE is" 1>&2
            exit 3 
        elif [ "\$CONTROLLED_ACCESS" = 'yes' ]; then
            echo "One or both of \$CDNA_FASTQ_FILE, \$BARCODES_FASTQ_FILE are not available at \$MANUAL_DOWNLOAD_FOLDER/ for this controlled access experiment" 1>&2
            exit 4   
        else
            confPart=''
            if [ -n "\$FASTQ_PROVIDER_CONFIG" ] && [ -e "\$FASTQ_PROVIDER_CONFIG" ]; then
                confPart=" -c \$FASTQ_PROVIDER_CONFIG"
            fi 

            # Stop fastq downloader from testing different methods -assume the control workflow has done that 
            export NOPROBE=1
        
            fetchFastq.sh -f ${WorkflowParamValidator.shellQuote(cdnaFastqURI)} -t "\$CDNA_FASTQ_FILE" -m "\$DOWNLOAD_METHOD" \$confPart
            
            # Allow for the first download also having produced the second output already

            if [ ! -e "\$BARCODES_FASTQ_FILE" ]; then
                fetchFastq.sh -f ${WorkflowParamValidator.shellQuote(barcodesFastqURI)} -t "\$BARCODES_FASTQ_FILE" -m "\$DOWNLOAD_METHOD" \$confPart
            fi
        fi
    """
}

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(WorkflowParamValidator.safeToken(row["${params.fields.run}"], params.fields.run), WorkflowParamValidator.safeToken(row["${params.fields.techrep}"], params.fields.techrep)) }
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

FINAL_FASTQS.into{
    FINAL_FASTQS_FOR_CONFIG
    FINAL_FASTQS_FOR_ALEVIN
}

// Derive Alevin barcodeconfig

process alevin_config {

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

    output:
        set val(runId), stdout into ALEVIN_CONFIG
    
    script:

        def barcodeConfig = ''

        if ( params.containsKey(protocol) ){

            canonicalProtocol = params.get(protocol)
            alevinType = canonicalProtocol.alevinType

            // Non-standard barcode config is supplied as a custom method

            if ( alevinType == 'custom' || "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
                barcodeConfig = "--barcodeLength ${barcodeLength} --umiLength ${umiLength} --end ${end}" 

            }else{
                barcodeConfig = "--$alevinType"
            }
            barcodeConfig = "-l ${canonicalProtocol.libType} $barcodeConfig" 
        }

        """
        if [ -z "$barcodeConfig" ]; then
            echo Input of $protocol results is misconfigured 1>&2
            exit 1
        fi

        # Also check barcode read lengths and return non-0 if they're not what they should be

        targetLen=\$(($umiLength + $barcodeLength))
        barcodesGood=0
        set +e
        while read -r l; do
            checkBarcodeRead.sh -r \$(readlink -f \$l) -b $barcodeLength -u $umiLength -n 1000000 1>&2
            if [ \$? -ne 0 ]; then
                barcodesGood=1
            fi
        done <<< "\$(ls barcodes*.fastq.gz)"
        set -e
        
        echo -n "$barcodeConfig"
        exit \$barcodesGood
        """
}

// Run Alevin per row
// Implement alevin_fry

process alevin {

    conda "${baseDir}/envs/alevin_fry.yml"
    
    cache 'deep'

    memory { 20.GB * task.attempt }
    cpus 12

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN.join(ALEVIN_CONFIG)
        file(transcriptToGene) from TRANSCRIPT_TO_GENE

    output:
        set val(runId), file("${runId}_ALEVIN_fry_quant") into ALEVIN_RESULTS
        set val(runId), file("${runId}_ALEVIN_fry_map/aux_info/meta_info.json") into ALEVIN_STATS

    script:

        canonicalProtocol = params.get(protocol)
        whitelist = canonicalProtocol.whitelist
        def protocolEnv = WorkflowParamValidator.shellQuote(params.protocol)
        def transcriptomeIndexEnv = WorkflowParamValidator.shellQuote(transcriptomeIndex)
        def transcriptToGeneEnv = WorkflowParamValidator.shellQuote(transcriptToGene)
        def runIdEnv = WorkflowParamValidator.shellQuote(runId)
        def whitelistEnv = WorkflowParamValidator.shellQuote(whitelist ?: '')
        def minMappingRateEnv = WorkflowParamValidator.shellQuote(params.minMappingRate)

        """
        PROTOCOL=${protocolEnv}
        TRANSCRIPTOME_INDEX=${transcriptomeIndexEnv}
        TRANSCRIPT_TO_GENE=${transcriptToGeneEnv}
        RUN_ID=${runIdEnv}
        WHITELIST=${whitelistEnv}
        MIN_MAPPING_RATE=${minMappingRateEnv}
        
        salmon alevin ${barcodeConfig} --sketch -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
            -i "\$TRANSCRIPTOME_INDEX" -p ${task.cpus} -o "\${RUN_ID}_ALEVIN_fry_map"

        if [ "\$PROTOCOL" = "10xv2" ]
        then
            alevin-fry generate-permit-list --input "\${RUN_ID}_ALEVIN_fry_map" -d fw --unfiltered-pl ${baseDir}/whitelist/737K-august-2016.txt --output-dir "\${RUN_ID}_ALEVIN_fry_quant_tmp" --min-reads 10
        elif [ "\$PROTOCOL" = "10xv3" ]
        then
            alevin-fry generate-permit-list --input "\${RUN_ID}_ALEVIN_fry_map" -d fw --unfiltered-pl "\$WHITELIST" --output-dir "\${RUN_ID}_ALEVIN_fry_quant_tmp" --min-reads 10
        elif [ "\$PROTOCOL" = "10x5prime" ]
        then
            alevin-fry generate-permit-list --input "\${RUN_ID}_ALEVIN_fry_map" -d rc --output-dir "\${RUN_ID}_ALEVIN_fry_quant_tmp" --force-cells 100000 --min-reads 10
        else
            alevin-fry generate-permit-list --input "\${RUN_ID}_ALEVIN_fry_map" -d fw --output-dir "\${RUN_ID}_ALEVIN_fry_quant_tmp" --force-cells 100000 --min-reads 10
        fi

        alevin-fry collate -i "\${RUN_ID}_ALEVIN_fry_quant_tmp" -r "\${RUN_ID}_ALEVIN_fry_map"
        alevin-fry quant -i "\${RUN_ID}_ALEVIN_fry_quant_tmp" -m "\$TRANSCRIPT_TO_GENE" -r cr-like-em -o "\${RUN_ID}_ALEVIN_fry_quant_tmp" --use-mtx

        TOTAL=\$(grep "num_processed" "\${RUN_ID}_ALEVIN_fry_map/aux_info/meta_info.json" |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')
        MAPPED=\$(grep "num_mapped" "\${RUN_ID}_ALEVIN_fry_map/aux_info/meta_info.json" |  awk '{split(\$0, array, ": "); print array[2]}'| sed 's/,//g')
        min_mapping=\$(echo "scale=2;((\$MAPPED * 100) / \$TOTAL)"|bc)

        if [ "\${min_mapping%.*}" -lt "\$MIN_MAPPING_RATE" ]; then
            echo "Minimum mapping rate (\$min_mapping) is less than the specified threshold of \$MIN_MAPPING_RATE" 1>&2
            exit 1 
        fi

        mv "\${RUN_ID}_ALEVIN_fry_quant_tmp" "\${RUN_ID}_ALEVIN_fry_quant"

        """

}

ALEVIN_RESULTS
    .into{
        ALEVIN_RESULTS_FOR_QC
        ALEVIN_RESULTS_FOR_PROCESSING
        ALEVIN_RESULTS_FOR_OUTPUT
    }

// Convert Alevin output to MTX. There will be one of these for every run, or
// technical replicate group of runs

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin_fry.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(runId), file(alevinResult) from ALEVIN_RESULTS_FOR_PROCESSING

    output:
        set val(runId), file("counts_mtx") into ALEVIN_MTX

    """
    RUN_ID=${WorkflowParamValidator.shellQuote(runId)}
    EXPERIMENT_TYPE=${WorkflowParamValidator.shellQuote(params.experimentType)}
    alevinFryMtxTo10x.py --cell_prefix "\${RUN_ID}-" "$alevinResult" counts_mtx "\$EXPERIMENT_TYPE"
    """ 
}

ALEVIN_MTX
    .into{
        ALEVIN_MTX_FOR_QC
        ALEVIN_MTX_FOR_EMPTYDROPS
        ALEVIN_MTX_FOR_OUTPUT
    }

// Make a diagnostic plot

ALEVIN_RESULTS_FOR_QC
    .join(ALEVIN_MTX_FOR_QC)
    .set{
        ALEVIN_QC_INPUTS
    }

process droplet_qc_plot{
    
    conda "${baseDir}/envs/droplet-barcode.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(runId), file(alevinResult), file(mtx) from ALEVIN_QC_INPUTS

    output:
        set val(runId), file("${runId}.png") into ALEVIN_QC_PLOTS

    """
    RUN_ID=${WorkflowParamValidator.shellQuote(runId)}
    dropletBarcodePlot.R --mtx-matrix counts_mtx/matrix.mtx --label "\$RUN_ID" --output-plot "\${RUN_ID}.png"
    """ 
}

// Remove empty droplets from Alevin results

process remove_empty_drops {
    
    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'ignore' }
    maxRetries 20
   
    input:
        set val(runId), file(countsMtx) from ALEVIN_MTX_FOR_EMPTYDROPS

    output:
        set val(runId), file('nonempty.rds') into NONEMPTY_RDS

    """
        EMPTY_DROPS_LOWER=${WorkflowParamValidator.shellQuote(params.emptyDrops.lower)}
        EMPTY_DROPS_NITERS=${WorkflowParamValidator.shellQuote(params.emptyDrops.nIters)}
        EMPTY_DROPS_FILTER_EMPTY=${WorkflowParamValidator.shellQuote(params.emptyDrops.filterEmpty)}
        EMPTY_DROPS_FILTER_FDR=${WorkflowParamValidator.shellQuote(params.emptyDrops.filterFdr)}
        MIN_CB_FREQ=${WorkflowParamValidator.shellQuote(params.minCbFreq)}
        dropletutils-read-10x-counts.R -s counts_mtx -c TRUE -o matrix.rds
        dropletutils-empty-drops.R -i matrix.rds --lower "\$EMPTY_DROPS_LOWER" --niters "\$EMPTY_DROPS_NITERS" --filter-empty "\$EMPTY_DROPS_FILTER_EMPTY" \
            --filter-fdr "\$EMPTY_DROPS_FILTER_FDR" --ignore "\$MIN_CB_FREQ" -o nonempty.rds -t nonempty.txt
    """
}

// Convert R matrix object with filtered cells back to .mtx

process rds_to_mtx{

    conda "${baseDir}/envs/dropletutils.yml"

    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
   
    input:
        set val(runId), file(rds) from NONEMPTY_RDS

    output:
        set val(runId), file("counts_mtx_nonempty") into NONEMPTY_MTX

    """ 
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))

        counts_sce <- readRDS('$rds')
        write10xCounts(assays(counts_sce)[[1]], path = 'counts_mtx_nonempty', barcodes = colData(counts_sce)\$Barcode, gene.id = rownames(counts_sce))
    """
}

// Compile raw results with raw and emptyDrops-filtered MTX

ALEVIN_RESULTS_FOR_OUTPUT
    .join(ALEVIN_MTX_FOR_OUTPUT)
    .join(NONEMPTY_MTX)
    .join(ALEVIN_QC_PLOTS)
    .join(ALEVIN_STATS)
    .set{ COMPILED_RESULTS }

process compile_results{

    publishDir "$resultsRoot/alevin", mode: 'copy', overwrite: true
    
    input:
        set val(runId), file('raw_alevin'), file(countsMtx), file(countsMtxNonempty), file(qcPlot), file(stats_file) from COMPILED_RESULTS

    output:
        set val(runId), file("$runId") into RESULTS_FOR_COUNTING

    """
        mkdir -p raw_alevin/alevin/mtx
        cp -P $countsMtx $countsMtxNonempty raw_alevin/alevin/mtx 
        mkdir -p raw_alevin/alevin/qc
        cp -P $qcPlot raw_alevin/alevin/qc
        cp -P raw_alevin $runId
        cp $stats_file $runId
    """
}

// Check the total number of runs we have 

RESULTS_FOR_COUNTING
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
