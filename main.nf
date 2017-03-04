#!/usr/bin/env nextflow

/* Preproccessing for variant calling
  Started October 2016

  @Authors
  Ryan Taylor <rwtaylor@stanford.edu>

  To run:
  $ nextflow run preprocessing.nf -c <nextflow.config> --sample <samples.tsv> --subsample [OPTIONAL] --validate [OPTIONAL]
  
  samples.tsv
  -----
  format: "sample library lane directory fastq1 fastq2"
  example: PT001 ng1 001 ../fastq-data PT001_TTCGGTG_L001_R1_001.fastq.gz PT001_TTCGGTG_L001_R2_001.fastq.gz

  NOTE: FreeBayes requires each read group to be unique. Or at least,
        that read groups are not shared between samples. Read group IDs
        will be created from "sample_lane_sampleGroup". So, it is
        important that each sample + lane + sampleGroup combination
        is unique to the read group.

  nextflow.config example
  -----
  process.executor = 'slurm'
  process.queue = 'dpetrov,hns,owners'

  params {
    sampleGroup  = 'cc1'
    genomeFasta  = '/scratch/PI/dpetrov/2016-tiger-wgs/reference/GCF_000464555.1_PanTig1.0_genomic.fna'
    genomeIndex  = '/scratch/PI/dpetrov/2016-tiger-wgs/reference/GCF_000464555.1_PanTig1.0_genomic.fna.fai'
    genomeDict   = '/scratch/PI/dpetrov/2016-tiger-wgs/reference/GCF_000464555.1_PanTig1.0_genomic.fna.dict'
    picardDir    = '/scratch/PI/dpetrov/local/opt'
    gatkDir      = '/scratch/PI/dpetrov/local/opt'
    params.popgenomics = '/scratch/PI/dpetrov/containers/singularity/popgenomics.img'
  }
  -----
*/

// Create fastq channel from samples.tsv

samples = file(params.samples)
fastqFiles = Channel
  .from(samples.readLines())
  .map {line ->
    lsplit      = line.split()
    sampleID    = lsplit[0]
    libID       = lsplit[1]
    laneID      = lsplit[2]
    fastqFile1  = file(lsplit[3]+ "/" + lsplit[4])
    fastqFile2  = file(lsplit[3]+ "/" + lsplit[5])
    [ sampleID, libID, laneID, fastqFile1, fastqFile2 ]
}

// Subsampling for testing
// To run:
// $ nextflow run preprocessing.nf -c <nextflow.config> --sample <samples.tsv> --subsample

if (params.subsample) {

  process MakeSmallTestData {

    publishDir "outputs/stages/subsampled-data"
    tag "${sampleID}-${libID}-${laneID}"

    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 5
    maxErrors '-1'

    input:
    set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles

    output:
    set sampleID, libID, laneID, file("${sampleID}_${libID}_${laneID}.small.1.fq"), file("${sampleID}_${libID}_${laneID}.small.2.fq") into fastqFilesSmall

    """
    seqtk sample -s100 ${fq1} 100000 > ${sampleID}_${libID}_${laneID}.small.1.fq
    seqtk sample -s100 ${fq2} 100000 > ${sampleID}_${libID}_${laneID}.small.2.fq

    """
  }

  fastqFilesSmall.into { fastqFiles_fastqc; fastqFiles_fastqvalidator; fastqFiles_trimgalore }

} else {

  fastqFiles.into { fastqFiles_fastqc; fastqFiles_fastqvalidator; fastqFiles_trimgalore }

}

process FastQC {
  publishDir "outputs/qc/fastqc", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { 1 * task.attempt }
  memory { 4.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles_fastqc

  output:
  file '*_fastqc.{zip,html}' into fastqc_results

  """
  /usr/local/bin/fastqc -q ${fq1} ${fq2}
  """
}

process Trim_galore {
  publishDir "outputs/stages/trimmed-fastqs"
  tag "${sampleID}-${libID}-${laneID}"

  cpus 2
  memory { 8.GB }
  time { 12.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles_trimgalore

  output:
  set sampleID, libID, laneID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs
  file '*trimming_report.txt' into trimgalore_results
  file '*_fastqc.{zip,html}' into fastqc_results_after_trimming

  """
  /usr/local/bin/trim_galore --paired --length 20 --gzip --fastqc ${fq1} ${fq2} 
  """
}

trimmedFastqs.into{ trimmedFastqs_mapping; trimmedFastqs_validate}

process Mapping {
  publishDir "outputs/stages/mapped-bams"
  tag "${sampleID}-${libID}-${laneID}"

  cpus 8
  memory { 32.GB * task.attempt }
  time { 24.h + (6.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file params.genomeFasta
  set sampleID, libID, laneID, file(fq1), file(fq2) from trimmedFastqs_mapping

  output:
  set sampleID, libID, laneID, file("*.bam") into mappedBams

  script:
  readGroupString="\"@RG\\tID:${sampleID}_${libID}_${laneID}\\tSM:${sampleID}\\tLB:${sampleID}_${libID}\\tPL:illumina\""

  """
  /usr/local/bin/bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} ${params.genomeFasta} ${fq1} ${fq2} | \
  /usr/local/bin/samtools view -hu - | /usr/local/bin/samtools sort -O bam - > ${sampleID}_${libID}_${laneID}.bam

  """
}

mappedBams.into { mappedBams_markduplicates; mappedBams_validate }

process MarkDuplicates_lane {
  publishDir "outputs/stages/marked-duplicates-bams-lane"
    tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 24.h + (6.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from mappedBams_markduplicates

  output:
  set sampleID, libID, laneID, file("*.md.bam"), file("*.md.bai") into markduplicatesBams
  set sampleID, libID, laneID, file("*.markduplicates.lane.txt") into markduplicatesResults

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar MarkDuplicates \
    TMP_DIR=picard_tmp \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.md.bam \
    M=${sampleID}_${libID}_${laneID}.markduplicates.lane.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true
  """
}

markduplicatesBams.into { markduplicatesBams; markduplicatesBams_index; markduplicatesBams_flagstat; markduplicatesBams_asMetrics; markduplicatesBams_isMetrics; markduplicatesBams_wgsMetrics ; markduplicatesBams_validate }

//  set sampleID, libID, laneID, file(bam), file("*.md.bai") into markduplicatesBams

markduplicatesBams_samples = markduplicatesBams.groupTuple(by:0)

// Combine bams with intervals

process MarkDuplicates_samples {
  publishDir "outputs/stages/marked-duplicates-bams-sample"
  tag "$sampleID"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 12.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bams), file(bais) from markduplicatesBams_samples
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  set sampleID, file("*.md.bam"), file("*.md.bai") into sampleBams
  set sampleID, file("*.markduplicates.samples.txt") into markduplicatesResults_samples

  script:
  input = bams.collect{"I=$it"}.join(' ')

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar MarkDuplicates \
    TMP_DIR=picard_tmp \
    $input \
    O=${sampleID}.md.bam \
    M=${sampleID}.markduplicates.samples.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.optical_duplicate_pixel_distance} \
    CREATE_INDEX=true
  """
}

sampleBams.into { sampleBams_flagstat; sampleBams_asMetrics; sampleBams_wgsMetrics; sampleBams_collectAlignment; sampleBams_preseq}

//
//
//
// ALIGNMENT METRICS ETC
// aligment is done, the rest is calculating metrics...
//
//
//
// LANE Metrics

process LaneCollectInsertSizeMetrics {
  publishDir "outputs/qc/lane-insert-size-metrics", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_isMetrics
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  file("${sampleID}_${libID}_${laneID}.txt") into lane_IsMetrics

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar CollectInsertSizeMetrics \
    TMP_DIR=picard_tmp \
    R=${gf} \
    I=${bam} \
    H=${sampleID}_${libID}_${laneID}.pdf \
    O=${sampleID}_${libID}_${laneID}.txt
  """
}

process LaneCollectAlignmentSummaryMetrics {
  publishDir "outputs/qc/lane-alignment-summary-metrics", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_asMetrics
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  file("${sampleID}_${libID}_${laneID}.txt") into lane_AsMetrics_lane

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar CollectAlignmentSummaryMetrics \
    TMP_DIR=picard_tmp \
    R=${gf} \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.txt
  """
}

process LaneCollectWgsMetrics {
  publishDir "outputs/qc/lane-wgs-metrics", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_wgsMetrics
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  file("${sampleID}_${libID}_${laneID}.txt") into lane_WgsMetrics

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar CollectWgsMetrics \
    TMP_DIR=picard_tmp \
    R=${gf} \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.txt
  """
}

process LaneFlagStat {
  publishDir "outputs/qc/lane-flagstat", mode: 'copy'
    tag "${sampleID}-${libID}-${laneID}"

  cpus { 2 * task.attempt }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'
  
  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_flagstat

  output:
  file("${sampleID}_${libID}_${laneID}.txt") into lane_flagstat

  """
  samtools flagstat ${sampleID}_${libID}_${laneID}.md.bam > ${sampleID}_${libID}_${laneID}.txt
  """
}


process SampleCollectAlignmentSummaryMetrics {
  publishDir "outputs/qc/sample-alignment-summary-metrics", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, file(bam) from sampleBams_asMetrics
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  file("${sampleID}.txt") into sample_AsMetrics

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar CollectAlignmentSummaryMetrics \
    TMP_DIR=picard_tmp \
    R=${gf} \
    I=${bam} \
    O=${sampleID}.txt
  """
}

process SampleCollectWgsMetrics {
  publishDir "outputs/qc/sample-wgs-metrics", mode: 'copy'
  tag "${sampleID}-${libID}-${laneID}"

  cpus { task.attempt == 1 ? 8: 16 }
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, file(bam) from sampleBams_wgsMetrics
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  file("${sampleID}.txt") into sample_WgsMetrics

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar CollectWgsMetrics \
    TMP_DIR=picard_tmp \
    R=${gf} \
    I=${bam} \
    O=${sampleID}.txt
  """
}

process SampleFlagStat {
  publishDir "outputs/qc/sample-flagstat", mode: 'copy'
  tag "$sampleID"

  cpus { 2 * task.attempt }
  time { 4.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, file(bam), file(bai) from sampleBams_flagstat

  output:
  file("${sampleID}.txt") into sample_flagstat

  """
  samtools flagstat ${sampleID}.md.real.bam > ${sampleID}.txt
  """
}

process SamplePreseq {
  publishDir "outputs/qc/sample-preseq", mode: 'copy'
  tag "$sampleID"
  
  cpus { 2 * task.attempt }
  time { 6.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'
  
  input:
  set sampleID, file(bam), file(bai) from sampleBams_preseq

  output:
  file("${sampleID}*") into sample_preseq

  """
  /usr/local/bin/preseq c_curve -o ${sampleID}.c_curve -bam ${bam}
  /usr/local/bin/preseq lc_extrap -o ${sampleID}.lc_extrap -bam ${bam}
  """
}

laneqc_files = lane_IsMetrics.mix(lane_AsMetrics_lane, lane_WgsMetrics, lane_flagstat).toList()
sampleqc_files = sample_AsMetrics.mix(sample_WgsMetrics, sample_flagstat, sample_preseq).toList()

process LaneMultiQC {

  publishDir "outputs", mode: 'copy'
  
  cpus { 2 * task.attempt }
  time { 6.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'
  
  input:
  file(laneqc_files) from laneqc_files

  output:
  file(laneqc) into laneqc_output

  """
  /usr/local/bin/multiqc -f -o laneqc .
  """
}

process SampleMultiQC {

  publishDir "outputs", mode: 'copy'
  
  cpus { 2 * task.attempt }
  time { 6.h + (2.h * task.attempt) }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
  maxRetries 5
  maxErrors '-1'
  
  input:
  file(sampleqc_files) from sampleqc_files

  output:
  file(sampleqc) into sampleqc_output

  """
  /usr/local/bin/multiqc -f -o sampleqc .
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

