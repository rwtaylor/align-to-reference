#!/usr/bin/env nextflow

/* Preproccessing for variant calling
  Started October 2016

  @Authors
  Ryan Taylor <rwtaylor@stanford.edu>

  To run:
  $ nextflow run preprocessing.nf -c <nextflow.config> --sample <samples.tsv> --subsample [OPTIONAL] --validate [OPTIONAL]
  
  samples.tsv
  -----
  format: "sample lane directory fastq1 fastq2"
  example: PT001 001 ../fastq-data PT001_TTCGGTG_L001_R1_001.fastq.gz PT001_TTCGGTG_L001_R2_001.fastq.gz

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

fastqFiles = fastqFiles.view()


// Subsampling for testing
// To run:
// $ nextflow run preprocessing.nf -c <nextflow.config> --sample <samples.tsv> --subsample

if (params.subsample) {

  process MakeSmallTestData {

    publishDir "outputs/subsampled-data"
    tag "${sampleID}-${libID}-${laneID}"

    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'

    input:
    set sampleID, libID, laneID, file(fqs), file(fq2) from fastqFiles

    output:
    set sampleID, libID, laneID, file("${sampleID}_${libID}_${laneID}.small.1.fq"), file("${sampleID}_${libID}_${laneID}.small.2.fq") into fastqFilesSmall

    """
    seqtk sample -s100 ${fq1} 10000 > ${sampleID}_${libID}_${laneID}.small.1.fq
    seqtk sample -s100 ${fq2} 10000 > ${sampleID}_${libID}_${laneID}.small.2.fq

    """
  }

  fastqFilesSmall = fastqFilesSmall.view()
  fastqFilesSmall.into { fastqFiles_fastqc; fastqFiles_fastqvalidator; fastqFiles_trimgalore }

} else {

  fastqFiles.into { fastqFiles_fastqc; fastqFiles_fastqvalidator; fastqFiles_trimgalore }

}

if (params.validate) {
 process FastQValidator {
   publishDir "outputs/validate-fq", mode: 'copy'
   tag "${sampleID}-${libID}-${laneID}"

   module "singularity"

   memory { 8.GB * task.attempt }
   time { 4.h * task.attempt }
   errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
   maxRetries 3
   maxErrors '-1'

   input:
   set idSample, idLane, file(fq1), file(fq2) from fastqFiles_fastqvalidator

   output:
   file '*.sanity.txt' into validated_fastqs

   """
   echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #" > ${idSample}_${idLane}.sanity.txt
   echo "FastQValidator and biawk Sanity Check" > ${idSample}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLane}.sanity.txt
   echo "FastQValidator results for ${fq1}" >> ${idSample}_${idLane}.sanity.txt
   singularity exec /scratch/PI/dpetrov/containers/singularity/bioinformatics.img fastQValidator --file ${fq1}  >> ${idSample}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLane}.sanity.txt
   echo "FastQValidator results for ${fq2}" >> ${idSample}_${idLane}.sanity.txt
   singularity exec /scratch/PI/dpetrov/containers/singularity/bioinformatics.img fastQValidator --file ${fq2} >> ${idSample}_${idLane}.sanity.txt

   echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #" >> ${idSample}_${idLane}.sanity.txt
   echo “Check if length of sequence equals quality score length”  >> ${idSample}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLane}.sanity.txt
   echo "Check results for ${fq1}" >> ${idSample}_${idLane}.sanity.txt
   zcat ${fq1} | bioawk -c fastx '{if (length(\$seq)!=length(\$qual)) print "Offending: " NR}' >> ${idSample}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLane}.sanity.txt
   echo "Check results for ${fq2}" >> ${idSample}_${idLane}.sanity.txt
   zcat ${fq2} | bioawk -c fastx '{if (length(\$seq)!=length(\$qual)) print "Offending: " NR}' >> ${idSample}_${idLane}.sanity.txt

   """
 }
}


process FastQC {
  publishDir "outputs/fastqc"
  tag "${sampleID}-${libID}-${laneID}"

  memory { 4.GB * task.attempt }
  time { 4.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles_fastqc

  output:
  file '*_fastqc.{zip,html}' into fastqc_results

  """
  source /scratch/PI/dpetrov/local/set_environment.sh
  fastqc -q ${fq1} ${fq2}
  """
}

process Trim_galore {
  publishDir "outputs/trimmed-fastqs"
  tag "${sampleID}-${libID}-${laneID}"

  cpus 4
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles_trimgalore

  output:
  set sampleID, libID, laneID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs
  file '*trimming_report.txt' into trimgalore_results

  """
  trim_galore --paired --length 20 --gzip ${fq1} ${fq2} 
  """
}

trimmedFastqs.into{ trimmedFastqs_mapping; trimmedFastqs_validate}

if (params.validate) {
 process FastQValidatorTrimmed {
   publishDir "outputs/validate-trimmed-fq", mode: 'copy'
   tag "${sampleID}-${libID}-${laneID}"

   module "singularity"

   memory { 8.GB * task.attempt }
   time { 4.h * task.attempt }
   errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
   maxRetries 3
   maxErrors '-1'

   input:
   set idSample, idLane, file(fq1), file(fq2) from trimmedFastqs_validate

   output:
   file '*.sanity.txt' into validated_trimmedfastqs

   """
   echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #" > ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "FastQValidator and biawk Sanity Check" > ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "FastQValidator results for ${fq1}" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   singularity exec /scratch/PI/dpetrov/containers/singularity/bioinformatics.img fastQValidator --file ${fq1}  >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "FastQValidator results for ${fq2}" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   singularity exec /scratch/PI/dpetrov/containers/singularity/bioinformatics.img fastQValidator --file ${fq2} >> ${idSample}_${idLib}_${idLane}.sanity.txt

   echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo “Check if length of sequence equals quality score length”  >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "Check results for ${fq1}" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   zcat ${fq1} | bioawk -c fastx '{if (length(\$seq)!=length(\$qual)) print "Offending: " NR}' >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   echo "Check results for ${fq2}" >> ${idSample}_${idLib}_${idLane}.sanity.txt
   zcat ${fq2} | bioawk -c fastx '{if (length(\$seq)!=length(\$qual)) print "Offending: " NR}' >> ${idSample}_${idLib}_${idLane}.sanity.txt

   """
 }
}

process Mapping {
  publishDir "outputs/mapped-bams"
  tag "${sampleID}-${libID}-${laneID}"

  cpus 8
  memory { 32.GB * task.attempt }
  time { 20.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  file params.genomeFasta
  set sampleID, libID, laneID, file(fq1), file(fq2) from trimmedFastqs_mapping

  output:
  set sampleID, libID, laneID, file("*.bam") into mappedBams

  script:
  readGroupString="\"@RG\\tID:${sampleID}_${libID}_${laneID}\\tSM:${sampleID}\\tLB:${sampleID}_${libID}\\tPL:illumina\""

  """
  bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} ${params.genomeFasta} ${fq1} ${fq2} | \
  samtools view -hu - | samtools sort -O bam - > ${sampleID}_${libID}_${laneID}.bam

  """
}

mappedBams.into { mappedBams_markduplicates; mappedBams_validate }

if (params.validate) {
 process ValidateBam {
   publishDir "outputs/validate-bams", mode: 'copy'
   tag "${sampleID}-${libID}-${laneID}"
   cpus { 8 * task.attempt }
   memory { 32.GB * task.attempt }
   time { 16.h * task.attempt }
   errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
   maxRetries 3
   maxErrors '-1'

   input:
   set idSample, idLane, file(bam) from mappedBams_validate

   output:
   file("*.bam.validate") into validated_mapped_bams

   """
   java -jar ${params.picardDir}/picard.jar ValidateSamFile \
     I=${sampleID}_${libID}_${laneID}.sorted.bam > ${sampleID}_${libID}_${laneID}.bam.validate

   """
 }
}

process MarkDuplicates {
  publishDir "outputs/marked-duplicates-bams"
    tag "${sampleID}-${libID}-${laneID}"

  cpus 8
  memory { 32.GB * task.attempt }
  time { 20.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from mappedBams_markduplicates

  output:
  set sampleID, libID, laneID, file("*.md.bam") into markduplicatesBams
  set file("*.markduplicates.txt") into markduplicatesResults

  """
  mkdir -p picard_tmp
  java -jar ${params.picardDir}/picard.jar MarkDuplicates \
    TMP_DIR=picard_tmp \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.md.bam \
    M=${sampleID}_${libID}_${laneID}.markduplicates.txt
  """
}

markduplicatesBams.into { markduplicatesBams; markduplicatesBams_index; markduplicatesBams_flagstat; markduplicatesBams_metrics; markduplicatesBams_wgsMetrics ; markduplicatesBams_validate }

process CollectMetrics {
  publishDir "outputs/picard-metrics"
  tag "${sampleID}-${libID}-${laneID}"

  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_metrics

  output:
  set sampleID, libID, laneID, file("${sampleID}.*") into metrics_md

  """
  mkdir -p picard_tmp
  java -jar ${params.picardDir}/picard.jar CollectMultipleMetrics \
    TMP_DIR=picard_tmp \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.metrics \
    PROGRAM=CollectQualityYieldMetrics
  """
}

process CollectWgsMetrics {
  publishDir "outputs/wgs-metrics"
  tag "${sampleID}-${libID}-${laneID}"

  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_wgsMetrics

  output:
  set sampleID, libID, laneID, file("${sampleID}.*") into wgsMetrics_md

  """
  mkdir -p picard_tmp
  java -jar ${params.picardDir}/picard.jar CollectWgsMetrics \
    TMP_DIR=picard_tmp \
    R=${params.genomeFasta} \
    I=${bam} \
    O=${sampleID}_${libID}_${laneID}.wgs_metrics
  """
}


process FlagStatMd {
  publishDir "outputs/flagstat-md"
    tag "${sampleID}-${libID}-${laneID}"

  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'
  
  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_flagstat

  output:
  file("*_md.flagstat.txt") into flagstat_md

  """
  samtools flagstat ${sampleID}_${libID}_${laneID}.md.bam > ${sampleID}_${libID}_${laneID}_md.flagstat.txt
  """
}

process IndexBams {
  tag "${sampleID}-${libID}-${laneID}"
  
  cpus 8
  memory { 32.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bam) from markduplicatesBams_index

  output:
  set sampleID, libID, laneID, file(bam), file("*.md.bam.bai") into markduplicatesBamsBais

  """
  samtools index ${bam}

  """
}

markduplicatesBamsBais = markduplicatesBamsBais.groupTuple(by:0)
markduplicatesBamsBais.into {mdbams_interval; mdbams_realignment}

process CreateIntervals {
  publishDir "outputs/bam-intervals"
  tag "$sampleID"

  cpus 8
  memory { 32.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(bams), file(bais) from mdbams_interval
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  set sampleID, file("*.intervals") into intervals

  script:
  input = bams.collect{"-I $it"}.join(' ')

  """
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkDir}/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $input \
  -R $gf \
  -nt ${task.cpus} \
  -o ${sampleID}.intervals

  """
}

// Combine bams with intervals
mdbams_realignment = mdbams_realignment.phase(intervals) {it -> it[0]}
mdbams_realignment = mdbams_realignment.map{a, b -> [a[0], a[1], a[2], a[3], a[4], b[1]]}

process Realign {
  publishDir "outputs/realigned-bams"
  tag "$sampleID"

  memory { 32.GB * task.attempt }
  time { 20.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, libIDs, laneIDs, file(bams), file(bais), file(intervals) from mdbams_realignment
  file gf from file(params.genomeFasta)
  file gi from file(params.genomeIndex)
  file gd from file(params.genomeDict)

  output:
  set sampleID, file("*.md.real.bam"), file("*.md.real.bai") into realignedBams

  script:
  input = bams.collect{"-I $it"}.join(' ')

  """
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkDir}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $input \
  -R $gf \
  -targetIntervals $intervals \
  -o '${sampleID}.md.real.bam'

  """
}

realignedBams.into { realignedBams_flagstat; realignedBams_metrics; realignedBams_wgsMetrics}

process FlagStatRealign {
  publishDir "outputs/flagstat-realigned"
  tag "$sampleID"

  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, file(bam), file(bai) from realignedBams_flagstat

  output:
  file("*_rln.flagstat.txt") into flagstat_realign

  """
  samtools flagstat ${sampleID}.md.real.bam > ${sampleID}_rln.flagstat.txt
  """
}

process CollectMetricsRealign {
  publishDir "outputs/picard-metrics-realigned"
  tag "$sampleID"
  
  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'
  
  input:
  set sampleID, file(bam), file(bai) from realignedBams_metrics

  output:
  set sampleID, libID, laneID, file("${sampleID}.*") into metrics_realign

  """
  mkdir -p picard_tmp
  java -jar ${params.picardDir}/picard.jar CollectMultipleMetrics \
        TMP_DIR=picard_tmp \
        I=${bam} \
        O=${sampleID}.metrics \
        PROGRAM=CollectQualityYieldMetrics
  """
}

process CollectWgsMetrics {
  publishDir "outputs/wgs-metrics-realigned"
  tag "$sampleID"
  
  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 3
  maxErrors '-1'

  input:
  set sampleID, file(bam), file(bai) from realignedBams_wgsMetrics

  output:
  set sampleID, file("${sampleID}.*") into realign_wgs_metrics

  """
  mkdir -p picard_tmp
  java -jar ${params.picardDir}/picard.jar CollectWgsMetrics \
        TMP_DIR=picard_tmp \
        R=${params.genomeFasta} \
        I=${bam} \
        O=${sampleID}.wgs_metrics
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

