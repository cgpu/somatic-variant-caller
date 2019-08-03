#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "            Somatic & Germline Variant Calling Pipeline             "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// fasta
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; ref_mutect2_tum_only_mode_channel ; ref_for_create_GenomicsDB_channel ; ref_create_somatic_PoN ; fasta_mutect; fasta_variant_eval  }
}
// fai
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
           .into { fai_mutect; fai_baserecalibrator; fai_haplotypecaller; ref_index_mutect2_tum_only_mode_channel ; ref_index_for_create_GenomicsDB_channel ; ref_index_create_somatic_PoN ; fai_variant_eval }
}

// dict
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_interval; dict_mutect; dict_baserecalibrator; dict_haplotypecaller; dict_variant_eval ; ref_dict_mutect2_tum_only_mode_channel ; ref_dict_for_create_GenomicsDB_channel ; ref_dict_create_somatic_PoN }
}

//dbsnp_gz
params.dbsnp_gz = params.genome ? params.genomes[ params.genome ].dbsnp_gz ?: false : false
if (params.dbsnp_gz) {
    Channel.fromPath(params.dbsnp_gz)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
           .set { dbsnp_gz}
}

//dbsnp_idx_gz
params.dbsnp_idx_gz = params.genome ? params.genomes[ params.genome ].dbsnp_idx_gz ?: false : false
if (params.dbsnp_idx_gz) {
    Channel.fromPath(params.dbsnp_idx_gz)
           .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
           .set { dbsnp_idx_gz}
}

//golden_indel_gz
params.golden_indel_gz = params.genome ? params.genomes[ params.genome ].golden_indel_gz ?: false : false
if (params.golden_indel_gz) {
    Channel.fromPath(params.golden_indel_gz)
           .ifEmpty { exit 1, "golden_indel_gz annotation file not found: ${params.golden_indel_gz}" }
           .set { golden_indel_gz }
}

//golden_indel_idx_gz
params.golden_indel_idx_gz = params.genome ? params.genomes[ params.genome ].golden_indel_idx_gz ?: false : false
if (params.golden_indel_idx_gz) {
    Channel.fromPath(params.golden_indel_idx_gz)
           .ifEmpty { exit 1, "golden_indel_idx_gz annotation file not found: ${params.golden_indel_idx_gz}" }
           .set { golden_indel_idx_gz }
}

//af_only_gnomad_vcf
params.af_only_gnomad_vcf = false
if (params.af_only_gnomad_vcf) {
    Channel.fromPath(params.af_only_gnomad_vcf)
           .ifEmpty { exit 1, "af_only_gnomad_vcf annotation file not found: ${params.af_only_gnomad_vcf}" }
           .into { af_only_gnomad_vcf_channel ; af_only_gnomad_vcf_channel_PoN }
}

//af_only_gnomad_vcf_idx
params.af_only_gnomad_vcf_idx = false
if (params.af_only_gnomad_vcf_idx) {
    Channel.fromPath(params.af_only_gnomad_vcf_idx)
           .ifEmpty { exit 1, "af_only_gnomad_vcf_idx annotation file not found: ${params.af_only_gnomad_vcf_idx}" }
           .into { af_only_gnomad_vcf_idx_channel ; af_only_gnomad_vcf_idx_channel_PoN }
}

//bwa_index_amb
params.bwa_index_amb = params.genome ? params.genomes[ params.genome ].bwa_index_amb ?: false : false
if (params.bwa_index_amb) {
    Channel.fromPath(params.bwa_index_amb)
           .ifEmpty { exit 1, "bwa_index_amb annotation file not found: ${params.bwa_index_amb}" }
           .set { bwa_index_amb }
}

//bwa_index_ann
params.bwa_index_ann = params.genome ? params.genomes[ params.genome ].bwa_index_ann ?: false : false
if (params.bwa_index_ann) {
    Channel.fromPath(params.bwa_index_ann)
           .ifEmpty { exit 1, "bwa_index_ann annotation file not found: ${params.bwa_index_ann}" }
           .set { bwa_index_ann }
}

//bwa_index_bwt
params.bwa_index_bwt = params.genome ? params.genomes[ params.genome ].bwa_index_bwt ?: false : false
if (params.bwa_index_bwt) {
    Channel.fromPath(params.bwa_index_bwt)
           .ifEmpty { exit 1, "bwa_index_bwt annotation file not found: ${params.bwa_index_bwt}" }
           .set { bwa_index_bwt }
}

//bwa_index_pac
params.bwa_index_pac = params.genome ? params.genomes[ params.genome ].bwa_index_pac ?: false : false
if (params.bwa_index_pac) {
    Channel.fromPath(params.bwa_index_pac)
           .ifEmpty { exit 1, "bwa_index_pac annotation file not found: ${params.bwa_index_pac}" }
           .set { bwa_index_pac }
}

//bwa_index_sa
params.bwa_index_sa = params.genome ? params.genomes[ params.genome ].bwa_index_sa ?: false : false
if (params.bwa_index_sa) {
    Channel.fromPath(params.bwa_index_sa)
           .ifEmpty { exit 1, "bwa_index_sa annotation file not found: ${params.bwa_index_sa}" }
           .set { bwa_index_sa }
}

// bed
if (params.bed) {
    Channel.fromPath(params.bed)
           .ifEmpty { exit 1, "BED file for to define regions not found: ${params.bed}" }
           .into { bed; bed_basename }

    bed_basename.map { file -> tuple(file.baseName, file) }.set{ bed_interval }
}

//intervals_list 
Channel.fromPath(params.interval_list_path, type: 'file')
       .into { intervals_haplotypecaller; intervals_mutect ;  interval_list_create_GenomicsDB_channel ; interval_list_mutect2_tum_only_mode_channel}

process gunzip_dbsnp {
    tag "$dbsnp_gz"

    input:
    file dbsnp_gz from dbsnp_gz
    file dbsnp_idx_gz from dbsnp_idx_gz

  	output:
  	file "*.vcf" into dbsnp, dbsnp_variantrecalibrator_snps, dbsnp_variantrecalibrator_indels
  	file "*.vcf.idx" into dbsnp_idx, dbsnp_idx_variantrecalibrator_snps, dbsnp_idx_variantrecalibrator_indels

    script:
    if ( "${dbsnp_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $dbsnp_gz
   	 gunzip -d --force $dbsnp_idx_gz
     """
   } else {
     """
     cp $dbsnp_gz dbsnp.vcf
     cp $dbsnp_idx_gz dbsnp.vcf.idx
     """
   }
}

process gunzip_golden_indel {
    tag "$golden_indel_gz"

    input:
    file golden_indel_gz from golden_indel_gz
    file golden_indel_idx_gz from golden_indel_idx_gz

    output:
    file "*.vcf" into golden_indel, golden_indel_variantrecalibrator_indels
    file "*.vcf.idx" into golden_indel_idx, golden_indel_idx_variantrecalibrator_indels

    script:
    if ( "${golden_indel_gz}".endsWith(".gz") ) {
        """
        gunzip -d --force $golden_indel_gz
        gunzip -d --force $golden_indel_idx_gz
        """
    } else {
        """
        cp $golden_indel_gz golden_indel.vcf
        cp $golden_indel_idx_gz golden_indel.vcf.idx
        """
    }
}

// TODO:
// Create junction for accepting BAM as input
Channel.fromPath(params.samples)
    .ifEmpty { exit 1, "samples file not found: ${params.samples}" }
    .splitCsv(sep: ',',  skip: 1 )
    .map{ shared_matched_pair_id, unique_subject_id, case_control_status, bam -> [shared_matched_pair_id, unique_subject_id, case_control_status, file(bam).baseName, file(bam)] }
    .into { samples; bams }

process BAM_sort {
    tag "$bam"
    container 'lifebitai/samtools:latest'

    input:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file(bam) from bams

    output:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file("${name}_mitoless.bam") into bam_sort, bam_sort_qc

    """
    samtools index $bam
    samtools view \
    -b $bam \
    chr21 chr22  > temp.bam && mv temp.bam ${name}.bam
    samtools sort -o ${name}_mitoless.bam ${name}.bam 
    rm ${name}.bam
    """
}

process RunBamQCmapped {
    tag "$bam"

    container 'maxulysse/sarek:latest'

    input:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file(bam) from bam_sort_qc

    output:
    file("${name}") into bamQCmappedReport

    when: !params.skip_multiqc

    script:
    // TODO: add --java-mem-size=${task.memory.toGiga()}G
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name} \
    -outformat HTML
    """
}

process MarkDuplicates {
    tag "$bam_sort"
    container 'broadinstitute/gatk:latest'

    input:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file(bam_sort) from bam_sort

    output:
    set val(name), file("${name}.bam"), file("${name}.bai"), val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status) into bam_markdup_baserecalibrator, bam_markdup_applybqsr
    file ("${name}.bam.metrics") into markDuplicatesReport

    """
    gatk MarkDuplicates \
    -I ${bam_sort} \
    --CREATE_INDEX true \
    -M ${name}.bam.metrics \
    -O ${name}.bam
    """
}

baserecalibrator_index = fasta_baserecalibrator.merge(fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = bam_markdup_baserecalibrator.combine(baserecalibrator_index)

process BaseRecalibrator {
    tag "$bam_markdup"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam_markdup), file(bai), val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), 
    file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

    output:
    set val(name), file("${name}_recal_data.table") into baserecalibrator_table
    file("*data.table") into baseRecalibratorReport

    """
    gatk BaseRecalibrator \
    -I $bam_markdup \
    --known-sites $dbsnp \
    --known-sites $golden_indel \
    -O ${name}_recal_data.table \
    -R $fasta
    """
}

applybqsr = baserecalibrator_table.join(bam_markdup_applybqsr)

process ApplyBQSR {
    tag "$baserecalibrator_table"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(baserecalibrator_table), file(bam), file(bai), val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status) from applybqsr

    output:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file("${name}_bqsr.bam"), file("${name}_bqsr.bai") into bam_for_qc, bam_haplotypecaller, bam_mutect

    script:
    """
    gatk ApplyBQSR -I $bam -bqsr $baserecalibrator_table -OBI -O ${name}_bqsr.bam
    """
}

process RunBamQCrecalibrated {
    tag "$bam"

    container 'maxulysse/sarek:latest'

    input:
    set val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file(bam), file(bai) from bam_for_qc

    output:
    file("${name}_recalibrated") into bamQCrecalibratedReport

    when: !params.skip_multiqc

    script:
    // TODO: add --java-mem-size=${task.memory.toGiga()}G \
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name}_recalibrated \
    -outformat HTML
    """
}

haplotypecaller_index = fasta_haplotypecaller.merge(fai_haplotypecaller, dict_haplotypecaller, bam_haplotypecaller)
haplotypecaller = intervals_haplotypecaller.combine(haplotypecaller_index)

process HaplotypeCaller {
    tag "$intervals_haplotypecaller"
    publishDir "${params.outdir}/GermlineVariantCalling", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    memory threadmem

    input:
    set file(intervals_haplotypecaller), file(fasta), file(fai), file(dict),
    val(shared_matched_pair_id), val(unique_subject_id), val(case_control_status), val(name), file(bam), file(bai) from haplotypecaller.collect()

    output:
    file("${name}.g.vcf") into haplotypecaller_gvcf
    file("${name}.g.vcf.idx") into index
    val(name) into name_mergevcfs

    when: !params.skip_haplotypecaller

    script:
    """
    gatk HaplotypeCaller \
    --java-options -Xmx${task.memory.toMega()}M \
    -R $fasta \
    -O ${name}.g.vcf \
    -I $bam \
    -ERC GVCF \
    -L $intervals_haplotypecaller
    """
}

bamsNormal = Channel.create()
bamsTumour = Channel.create()

bam_mutect.choice( bamsTumour, bamsNormal ) { it[2] =~ 1 ? 0 : 1 }
bamsNormal.into {bamsNormal_PoN ; bamsNormal_mutect}
combined_bam = bamsNormal_mutect.combine(bamsTumour, by: 0)

bamsNormal_PoN.into {bamsNormal_PoN_bam_ ; bamsNormal_PoN_bai_ }
bamsNormal_PoN_bam = bamsNormal_PoN_bam_.map { shared_matched_pair_id, unique_subject_id, case_control_status, name, bam, bai -> [bam]}
bamsNormal_PoN_bai = bamsNormal_PoN_bai_.map { shared_matched_pair_id, unique_subject_id, case_control_status, name, bam, bai -> [bai]}

ref_mutect = fasta_mutect.merge(fai_mutect, dict_mutect)
variant_calling = combined_bam.combine(ref_mutect)
variant_calling_intervals = intervals_mutect.combine(variant_calling)
variant_calling_intervals.into{ mutect; names_for_vcf2maf}


process run_mutect2_tumor_only_mode {

    tag "${normal_bam.simpleName.minus('_Normal').minus('_bqsr')}"
    publishDir "${params.outdir}/MutectTumorOnlyMode", mode: 'copy'
    container "broadinstitute/gatk:latest"

    input:
    file(normal_bam) from bamsNormal_PoN_bam.collect()
    file(normal_bai) from bamsNormal_PoN_bai.collect()
    each file(ref) from ref_mutect2_tum_only_mode_channel
    each file(ref_index) from ref_index_mutect2_tum_only_mode_channel
    each file(ref_dict) from ref_dict_mutect2_tum_only_mode_channel
    each file(intervals) from interval_list_mutect2_tum_only_mode_channel

    output:
    file('*.vcf.gz') into vcf_for_create_GenomicsDB_channel
    file('*.vcf.gz.tbi') into vcf_tbi_for_create_GenomicsDB_channel

    script:
    """

    gatk Mutect2 \
    -R ${ref} \
    -I ${normal_bam} -normal ${normal_bam.simpleName.minus('_Normal').minus('_bqsr')} \
    --max-mnp-distance 0 \
    -O ${normal_bam.baseName}.vcf.gz \
    -L $intervals \
    --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
    """
}

process create_GenomicsDB {

    tag "all_the_vcfs"
    publishDir "${params.outdir}/GenomicsDBImport", mode: 'copy'
    container "broadinstitute/gatk:latest"

    input:
    file("*") from vcf_for_create_GenomicsDB_channel.collect()
    file("*") from vcf_tbi_for_create_GenomicsDB_channel.collect()
    file(ref) from ref_for_create_GenomicsDB_channel
    file(ref_index) from ref_index_for_create_GenomicsDB_channel
    file(ref_dict) from ref_dict_for_create_GenomicsDB_channel
    file(intervals) from interval_list_create_GenomicsDB_channel

    output:
    file("pon_db") into pon_db_for_create_somatic_PoN

    shell:
    '''
    echo -n "gatk GenomicsDBImport -R !{ref} --genomicsdb-workspace-path pon_db " > create_GenomicsDB.sh
    for vcf in $(ls *.vcf.gz); do
    echo -n "-V $vcf " >> create_GenomicsDB.sh
    done
    echo -n "-L !{intervals}" --merge-input-intervals --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' >> create_GenomicsDB.sh
    chmod ugo+xr create_GenomicsDB.sh
    bash create_GenomicsDB.sh
    chmod -R ugo+xrw pon_db
    '''
}

process create_somatic_PoN {
    
    tag "$af_only_gnomad_vcf"
    publishDir "${params.outdir}/CreateSomaticPanelOfNormals", mode: 'copy'
    container "broadinstitute/gatk:latest"

    input:
    file(pon_db) from pon_db_for_create_somatic_PoN
    file(ref) from ref_create_somatic_PoN
    file(ref_index) from ref_index_create_somatic_PoN
    file(ref_dict) from ref_dict_create_somatic_PoN
    file(af_only_gnomad_vcf) from af_only_gnomad_vcf_channel_PoN
    file(af_only_gnomad_vcf_idx) from af_only_gnomad_vcf_idx_channel_PoN

    output:
    set file("pon.vcf.gz"), file("pon.vcf.gz.tbi") into create_somatic_PoN_results_channel
    
    script:
    """
    gatk CreateSomaticPanelOfNormals \
    -R $ref \
    --germline-resource $af_only_gnomad_vcf \
    -V gendb://$pon_db \
    -O pon.vcf.gz  
    """
}

process Mutect2 {

    tag "${tumourSampleId}_vs_${sampleId}.vcf"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/Somatic", mode: 'copy'

    input:
    set file(intervals_mutect), val(patientId), val(sampleId), val(status), val(name), file(bam), file(bai),
    val(tumourSampleId), val(tumourStatus), val(tumourName), file(tumourBam), file(tumourBai),
    file(fasta), file(fai), file(dict) from mutect.collect()
    file(af_only_gnomad_vcf) from af_only_gnomad_vcf_channel
    file(af_only_gnomad_vcf_idx) from af_only_gnomad_vcf_idx_channel
    set file(pon_vcf_gz), file(pon_vcf_gz_tbi) from create_somatic_PoN_results_channel

    output:
    set val("${tumourSampleId}_vs_${sampleId}"), file("${tumourSampleId}_vs_${sampleId}.vcf"),  file("${tumourSampleId}_vs_${sampleId}.vcf.idx"), file("${tumourSampleId}_vs_${sampleId}.vcf.stats") into vcf_variant_eval, vcf_for_vcf2maf, vcf_for_filter_mutect_calls

    script:
    """
    tumourName_bash=`echo ${tumourName}`
    name_bash=`echo ${name}`

    tumourName_trimmed=`echo \${tumourName_bash%_*}`
    name_trimmed=`echo \${name_bash%_*}`

    gatk Mutect2 \
    -R ${fasta} \
    -I ${tumourBam}  -tumor \${tumourName_trimmed} \
    -I ${bam} -normal \${name_trimmed} \
    -O ${tumourSampleId}_vs_${sampleId}.vcf \
    -L $intervals_mutect \
    --panel-of-normals  $pon_vcf_gz \
    --germline-resource $af_only_gnomad_vcf \
    --interval-padding 100 
    #gatk --java-options "-Xmx\${task.memory.toGiga()}g" \
    """
}

vcf2maf = vcf_for_vcf2maf.combine(names_for_vcf2maf)


process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    container 'ewels/multiqc:v1.7'

    when:
    !params.skip_multiqc

    input:
    file (bam_metrics) from markDuplicatesReport.collect().ifEmpty([])
    file (bamQC) from bamQCmappedReport.collect().ifEmpty([])
    file (bamQCrecalibrated) from bamQCrecalibratedReport.collect().ifEmpty([])
    file (baseRecalibrator) from baseRecalibratorReport.collect().ifEmpty([])
    
    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . -m picard -m qualimap -m gatk
    """
}