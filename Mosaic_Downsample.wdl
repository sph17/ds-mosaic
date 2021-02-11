version 1.0

import "Structs.wdl"
import "downsampling_mosaic_part1.wdl" as ds1
import "downsampling_mosaic_part2.wdl" as ds2

workflow downSampling_mosaic {

  #################################################################################
  ####        Required basic arguments for downsampling pipeline                  #
  #################################################################################
  
  input {
    File reference_fasta
    String downsample_docker
    File bam_or_cram_file_1
    File bam_or_cram_file_2

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    File ref_dict

    String sample_ID
    File intervals_genome

    #default start depth is 30x
    Float start_depth_1 = 30
    Float? start_depth_2
    Int? seed_override  

    File? gatk4_jar_override
    String gatk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_bam_to_fastq 
    RuntimeAttr? runtime_attr_random_sample
    RuntimeAttr? runtime_attr_realign
    RuntimeAttr? runtime_attr_add_read_group
    RuntimeAttr? runtime_attr_mark_duplicates
    RuntimeAttr? runtime_attr_sort_index
    RuntimeAttr? runtime_attr_count_coverage
    RuntimeAttr? runtime_attr_collect_counts
    RuntimeAttr? runtime_attr_combine_random_sort_1  
    RuntimeAttr? runtime_attr_combine_random_sort_2

    #Execution defaults and overrides
    Boolean run_downsample_custom = true
    Boolean process_downsample_custom = true

    Float? final_depth_custom_1
    Float? final_depth_custom_2

  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
  
  #################################################################################
  ####        Calls part 1 to convert cram/bam file to paired fastq               #
  #################################################################################
  if (run_downsample_custom){
    call ds1.downSampling_01 as downSampling_1{
      input :
        bam_or_cram_file = bam_or_cram_file_1,
        reference_fasta = reference_fasta,
        downsample_docker = downsample_docker,
        seed_override=seed_override,
        start_depth = start_depth_1,
        final_depth = select_first([final_depth_custom_1,]),
        runtime_attr_cram_to_bam = runtime_attr_cram_to_bam,
        runtime_attr_bam_to_fastq = runtime_attr_bam_to_fastq,
        runtime_attr_random_sample = runtime_attr_random_sample
    }
  }
  

  #################################################################################
  ####        Calls part 1 to convert cram/bam file to paired fastq               #
  #################################################################################
  if (run_downsample_custom){
    call ds1.downSampling_01 as downSampling_2{
      input :
        bam_or_cram_file = bam_or_cram_file_2,
        reference_fasta = reference_fasta,
        downsample_docker = downsample_docker,
        seed_override = seed_override,
        start_depth = select_first([start_depth_2,start_depth_1]),
        final_depth = select_first([final_depth_custom_2,]),
        runtime_attr_cram_to_bam = runtime_attr_cram_to_bam,
        runtime_attr_bam_to_fastq = runtime_attr_bam_to_fastq,
        runtime_attr_random_sample = runtime_attr_random_sample 
    }
  }


  #################################################################################
  ####        Downsamples samples to custom read depths                           #
  #################################################################################

  if (process_downsample_custom && run_downsample_custom) {
    call ds2.downSampling_02 {
      input :
        fastq_1_file_1 = select_first([downSampling_1.fastq_1,]), #.1.fastq file for first sample fastq file
        fastq_1_file_2 = select_first([downSampling_1.fastq_2,]), #.2.fastq file for first sample fastq file
        fastq_2_file_1 = select_first([downSampling_2.fastq_1,]), #.1.fastq file for second sample fastq file
        fastq_2_file_2 = select_first([downSampling_2.fastq_2,]), #.1.fastq file for second sample fastq file
        downsample_docker = downsample_docker,
        reference_fasta = reference_fasta,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        original_cram_or_bam_file_read_groups = select_first([downSampling_1.read_groups,]), #puts back readgroup from first sample file
        intervals_genome = intervals_genome,
        sample_ID = sample_ID, #the name you want to call your new in silico mix sample
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        runtime_attr_realign = runtime_attr_realign,
        runtime_attr_add_read_group = runtime_attr_add_read_group,
        runtime_attr_mark_duplicates = runtime_attr_mark_duplicates,
        runtime_attr_sort_index = runtime_attr_sort_index,
        runtime_attr_count_coverage = runtime_attr_count_coverage,
        runtime_attr_collect_counts = runtime_attr_collect_counts
    }
  }

  output {
    #sample 1/A
    File? fastq_1A = downSampling_1.fastq_1
    File? fastq_2A = downSampling_1.fastq_2
    File? read_groups_A = downSampling_1.read_groups
    #sample 2/B
    File? fastq_1B = downSampling_2.fastq_1
    File? fastq_2B = downSampling_2.fastq_2
    File? read_groups_B = downSampling_2.read_groups

    #sample in-silico mix
    File? markdup_metrics_custom = downSampling_02.markdup_metrics
    File? sorted_cram_custom = downSampling_02.sorted_cram
    File? crai_file_custom = downSampling_02.crai_file
    File? wgs_coverage_metrics_custom = downSampling_02.wgs_coverage_metrics
    File? read_counts_custom = downSampling_02.read_counts
  }

}



