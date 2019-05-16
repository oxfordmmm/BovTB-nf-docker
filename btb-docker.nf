#!/usr/bin/env nextflow

params.input_dir = "/data/tests/input/"
params.reads = "*_{1,2}.fastq.gz"
data_path = params.input_dir + params.reads
params.output_dir = "/data/tests/output"
params.lowmem = ""
lowmem = Channel.value("${params.lowmem}")

/* location of reference information */
ref = file(params.ref)
refgbk = file(params.refgbk)
stage1pat = file(params.stage1pat)
stage2pat = file(params.stage2pat)
adapters = file(params.adapters)
kraken2db = file(params.kraken2db)

/*	Collect pairs of fastq files and infer sample names */
Channel
    .fromFilePairs( data_path, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 
	read_pairs.into { read_pairs; raw_reads }


/* remove duplicates from raw data */
process Deduplicate {
    tag {pair_id}
	
	maxForks 2

	publishDir "${params.output_dir}/deduplicate", mode: "copy"

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz") from read_pairs

	output:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into dedup_read_pairs
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into uniq_reads

	"""
	gunzip -c ${pair_id}_*_R1_*.fastq.gz > ${pair_id}_R1.fastq 
	gunzip -c ${pair_id}_*_R2_*.fastq.gz > ${pair_id}_R2.fastq
	echo '${pair_id}_R1.fastq\n${pair_id}_R2.fastq' > fqin.lst
	${FASTUNIQ}/source/fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
	rm ${pair_id}_R1.fastq
	rm ${pair_id}_R2.fastq
	"""
}	

/* trim adapters and low quality bases from fastq data */
/* ILLUMINACLIP:/home/richard/ReferenceSequences/adapter.fasta:2:30:10 \ */
process Trim {
	tag {pair_id}

	maxForks 2

	publishDir "${params.output_dir}/trim", mode: "copy"

	input:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") from dedup_read_pairs

	output:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs2
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_reads
	
	"""
	java -jar ${TRIM}/trimmomatic-0.38.jar PE -threads 3 -phred33 \
	${pair_id}_uniq_R1.fastq ${pair_id}_uniq_R2.fastq  \
	${pair_id}_trim_R1.fastq ${pair_id}_fail1.fastq \
	${pair_id}_trim_R2.fastq ${pair_id}_fail2.fastq \
	ILLUMINACLIP:${adapters}:2:30:10 \
	SLIDINGWINDOW:10:20 MINLEN:36 

	rm ${pair_id}_fail1.fastq
	rm ${pair_id}_fail2.fastq
	"""
}


// Building index for ref fasta
process BWA_Index {
    tag {ref}

    publishDir "${params.output_dir}/indexed_ref", mode: "copy"

    input:
        file ref

    output:
        file "*" into bwa_index

    script:
    """
    $BWA/bwa index ${ref}
    """
}

/* map to reference sequence */
process Map2Ref {
    tag {pair_id}

	maxForks 1

	publishDir "${params.output_dir}/map2ref", mode: "copy"

	input:
	file ref
	file index from bwa_index
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from trim_read_pairs

	output:
	set pair_id, file("${pair_id}.mapped.sorted.bam") into mapped_bam
	set pair_id, file("${pair_id}.mapped.sorted.bam") into bam4stats

	"""
	${BWA}/bwa mem -T10 -M -t2 ${ref}  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq \
	| ${SAMTOOLS}/samtools view -@2 -ShuF 2308 - \
	| ${SAMTOOLS}/samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam
	"""
}


/* Variant calling */
process VarCall {
    tag {pair_id}

	maxForks 4

	publishDir "${params.output_dir}/varcall", mode: "copy"

	input:
	file ref
	set pair_id, file("${pair_id}.mapped.sorted.bam") from mapped_bam

	output:
	set pair_id, file("${pair_id}.pileup.vcf.gz") into vcf
	set pair_id, file("${pair_id}.pileup.vcf.gz") into vcf2
	set pair_id, file("${pair_id}.pileup.vcf.gz") into vcf3

	"""
	${SAMTOOLS}/samtools index ${pair_id}.mapped.sorted.bam	
	${BCFTOOLS}/bcftools mpileup -q 60 -Ou -f $ref ${pair_id}.mapped.sorted.bam \
	| ${BCFTOOLS}/bcftools call --ploidy 1 -cf GQ - -Oz -o  ${pair_id}.pileup.vcf.gz
	"""
}

/* Consensus calling */
process VCF2Consensus {
	tag {pair_id}
	
	publishDir "${params.output_dir}/consensus", mode: 'copy', pattern: '*_consensus.fas'
	publishDir "${params.output_dir}/bcf", mode: 'copy', pattern: '*.norm-flt.bcf'

	maxForks 2

	input:
	set pair_id, file("${pair_id}.pileup.vcf.gz") from vcf2

	output:
	set pair_id, file("${pair_id}_consensus.fas"), file("${pair_id}.norm-flt.bcf") into consensus

	"""
	${BCFTOOLS}/bcftools index ${pair_id}.pileup.vcf.gz
	${BCFTOOLS}/bcftools norm -f $ref ${pair_id}.pileup.vcf.gz -Ob | ${BCFTOOLS}/bcftools filter --IndelGap 5 - -Ob -o ${pair_id}.norm-flt.bcf
	${BCFTOOLS}/bcftools index ${pair_id}.norm-flt.bcf
	${BCFTOOLS}/bcftools consensus -f $ref ${pair_id}.norm-flt.bcf | sed '/^>/ s/.*/>${pair_id}/' - > ${pair_id}_consensus.fas
	"""
}

//	Combine data for generating per sample statistics

raw_reads
	.join(uniq_reads)
	.set { raw_uniq }

trim_reads
	.join(bam4stats)
	.set { trim_bam }

raw_uniq
	.join(trim_bam)
	.set { input4stats }

/* Mapping Statistics*/
process ReadStats{
    tag {pair_id}

	maxForks 2

	publishDir "${params.output_dir}/readstats", mode: "copy"

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz"), file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq"), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq"), file("${pair_id}.mapped.sorted.bam") from input4stats

	output:
	set pair_id, file("${pair_id}_stats.csv") into stats
	set pair_id, file('outcome.txt') into Outcome

	shell:
	'''
	raw_R1=$(zgrep -c "^+$" !{pair_id}_*_R1_*.fastq.gz)
	uniq_R1=$(grep -c "^+$" !{pair_id}_uniq_R1.fastq)
	trim_R1=$(grep -c "^+$" !{pair_id}_trim_R1.fastq)
	num_map=$(${SAMTOOLS}/samtools view -c !{pair_id}.mapped.sorted.bam)
	avg_depth=$(${SAMTOOLS}/samtools depth  !{pair_id}.mapped.sorted.bam  |  awk '{sum+=$3} END { print sum/NR}')

	num_raw=$(($raw_R1*2))
	num_uniq=$(($uniq_R1*2))
	num_trim=$(($trim_R1*2))
	pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
	pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_raw)" |bc)
	pc_mapped=$(echo "scale=2; ($num_map*100/$num_trim)" |bc)

	mindepth=10
	minpc=60
	minreads=600000
	
	if [ ${avg_depth%%.*} -ge $mindepth ] && [ ${pc_mapped%%.*} -gt $minpc ]; then flag="Pass"
		elif [ ${avg_depth%%.*} -lt $mindepth ] && [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="Comtaminated"
		elif [ ${avg_depth%%.*} -lt $mindepth ] && [ $num_trim -lt $minreads ]; then flag="InsufficientData"
#		elif [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="q_OtherMycobact"
		else flag="CheckRequired"
	fi
 
	echo "Sample,NumRawReads,NumDedupReads,%afterDedup,NumTrimReads,%afterTrim,NumMappedReads,%Mapped,MeanCov,Outcome" > !{pair_id}_stats.csv
	echo "!{pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$num_map","$pc_mapped","$avg_depth","$flag"" >> !{pair_id}_stats.csv
	echo "$flag" > outcome.txt
	'''	
}

/* SNP filtering and annotation */
process SNPfiltAnnot{
    tag {pair_id}

	maxForks 4

	publishDir "${params.output_dir}/snpfiltannot", mode: "copy"

	input:
	file ref
	file refgbk
	set pair_id, file("${pair_id}.pileup.vcf.gz") from vcf3

	output:
	set pair_id, file("${pair_id}.pileup_SN.csv"), file("${pair_id}.pileup_DUO.csv"), file("${pair_id}.pileup_INDEL.csv") into VarTables
	set pair_id, file("${pair_id}.pileup_SN_Annotation.csv") into VarAnnotation

	"""
	${BCFTOOLS}/bcftools view -O v ${pair_id}.pileup.vcf.gz \
	| python ${SCRIPT_PATH}/snpsFilter.py - ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp}
	mv _DUO.csv ${pair_id}.pileup_DUO.csv
	mv _INDEL.csv ${pair_id}.pileup_INDEL.csv
	mv _SN.csv ${pair_id}.pileup_SN.csv
	python ${SCRIPT_PATH}/annotateSNPs.py ${pair_id}.pileup_SN.csv ${refgbk} ${ref}
	"""
}

/*	Combine data for assign cluster for each sample*/
vcf
	.join(stats)
	.set { input4Assign }

/* Assigns cluster by matching patterns of cluster specific SNPs. Also suggests inferred historical genotype */
process AssignClusterCSS{
    tag {pair_id}

	maxForks 2

	publishDir "${params.output_dir}/assignclustercss", mode: "copy"

	input:
	file ref
	set pair_id, file("${pair_id}.pileup.vcf.gz"), file("${pair_id}_stats.csv") from input4Assign

	output:
	file("${pair_id}_stage1.csv") into AssignCluster

	"""
	gunzip -c ${pair_id}.pileup.vcf.gz > ${pair_id}.pileup.vcf
	python ${SCRIPT_PATH}/Stage1-test.py ${pair_id}_stats.csv ${stage1pat} ${ref} test 1 ${min_mean_cov} ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${min_qual_nonsnp} ${pair_id}.pileup.vcf
	mv _stage1.csv ${pair_id}_stage1.csv
	"""
}

Outcome
	.join(trim_read_pairs2)
	.set {IDdata}

/* Identify any non-M.bovis samples using kraken */
process IDnonbovis{

	tag {pair_id}

	publishDir "${params.output_dir}/kraken2", mode: 'copy', pattern: '*.tab'

	maxForks 1

	input:
	set pair_id, file('outcome.txt'), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from IDdata
	val lowmem from lowmem

	output:
	set pair_id, file("${pair_id}_kraken2.tab") into IDnonbovis

	"""
	${KRAKEN2}/kraken2 --threads 2 --quick $lowmem --db $kraken2db --output - --report ${pair_id}_kraken2.tab --paired ${pair_id}_trim_R1.fastq  ${pair_id}_trim_R2.fastq 
	"""
}


/* Combine all cluster assignment data into a single results file */
AssignCluster
	.collectFile( name: 'AssignedWGSCluster.csv', sort: true, storeDir: "${params.output_dir}/AssignClusterCSS", keepHeader: true )


workflow.onComplete {
	log.info "Completed sucessfully:	$workflow.success"		
	log.info "Nextflow Version:	$workflow.nextflow.version"
	log.info "Duration:		$workflow.duration"
	log.info "Output Directory:	$params.output_dir"
}






