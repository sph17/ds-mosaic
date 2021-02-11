version 1.0

import "Structs.wdl"

workflow downSampling_01 {

  #################################################################################
  ####        Required basic arguments for downsampling part 1                    #
  #################################################################################
    
  input {
    File bam_or_cram_file
    File reference_fasta
    String downsample_docker
    Float start_depth
    Float final_depth
    Int? seed_override

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fai
    File ref_dict

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_bam_to_fastq
    RuntimeAttr? runtime_attr_random_sample    
  }

  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
    bam_file: ".bam file to search for SVs. bams are preferable, crams will be converted to bams."
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    start_depth: "float/integer for initial depth of sequencing."
    final_depth: "float/integer for final desired depth."
    seed: "optional seed input for random sampling."
  }

  meta {
    author: "Stephanie Hao"
    email: "shao@broadinstitute.org"
  }
    
  #################################################################################
  ####        Convert cram to bam                                                 #
  #################################################################################
  
  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  
  # Convert to BAM if we have a CRAM
  if (!is_bam_) {
    call cramToBam {
      input :
        cram_file = bam_or_cram_file,
        reference_fasta = reference_fasta,
        downsample_docker = downsample_docker,
        runtime_attr_override = runtime_attr_cram_to_bam
    }
  }


  #################################################################################
  ####        Convert bam to fq1 & fq2                                            #
  #################################################################################
  call bamToFq { 
    input :
      bam_file = select_first([cramToBam.bam_file, bam_or_cram_file]),
      downsample_docker = downsample_docker,
      runtime_attr_override = runtime_attr_bam_to_fastq
  }

  #################################################################################
  ####        Calculate reads and random samples appropriate # of reads           #
  #################################################################################
  call countAndRandomSample {
    input : 
      fastq_file_1 = bamToFq.fastq_file_1,
      fastq_file_2 = bamToFq.fastq_file_2,
      downsample_docker = downsample_docker,
      start_depth = start_depth,
      final_depth = final_depth,
      seed = select_first([seed_override, 20937]),
      runtime_attr_override = runtime_attr_random_sample
    }


  output {
    File fastq_1 = bamToFq.fastq_file_1
    File fastq_2 = bamToFq.fastq_file_2
    File read_groups = bamToFq.read_groups_file
    File downsample_file_1 = countAndRandomSample.downsample_file_1
    File downsample_file_2 = countAndRandomSample.downsample_file_2
  }

}


task cramToBam {
    
  input {
    File cram_file
    File reference_fasta
    String downsample_docker
    RuntimeAttr? runtime_attr_override
  }

  File reference_index_file = reference_fasta + ".fai"
    
  String bam_file_name = basename(cram_file, ".cram") + ".bam"

  Int num_cpu = 1
  Float mem_size_gb = num_cpu * 4.0
  
  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 20.0
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(bam_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bam_file = bam_file_name
  }

  command {
    
    set -euo pipefail

    #converts cram files to bam files
    samtools view \
            -b \
            -h \
            -@ ~{num_cpu} \
            -T "~{reference_fasta}" \
            -o "~{bam_file_name}" \
            "~{cram_file}"
    }

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task bamToFq {
  
  input {
    File bam_file
    String downsample_docker
    RuntimeAttr? runtime_attr_override
  }

  String fastq_file_1_name = basename(bam_file, ".bam") + "_1.fastq"
  String fastq_file_2_name = basename(bam_file, ".bam") + "_2.fastq"
  String read_groups_name = basename(bam_file, ".bam") + "_read_groups.txt"

  Int num_cpu = 5
  Int mem_size_gb = 16
  Int vm_disk_size = 300

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File fastq_file_1 = fastq_file_1_name
    File fastq_file_2 = fastq_file_2_name
    File read_groups_file = read_groups_name
  }

  command <<<
    set -euo pipefail

    #extracts read groups for part 2 of pipeline
    samtools view -H ~{bam_file} > ~{read_groups_name}
    

    #converts bam file to paired fastq files
    java -Xmx16G -jar /opt/conda/share/picard-2.23.8-0/picard.jar SamToFastq \
            I=~{bam_file} \
            FASTQ=~{fastq_file_1_name} \
            SECOND_END_FASTQ=~{fastq_file_2_name}
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task countAndRandomSample {
  
  input {
    File fastq_file_1
    File fastq_file_2
    Float start_depth
    Float final_depth
    String downsample_docker
    Int seed
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = 1
  Int mem_size_gb = 6
  Int vm_disk_size = 400

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }  

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String fastq_downsample_1_name = basename(fastq_file_1, "_1.fastq") + "_downsample_~{final_depth}x"
  String fastq_downsample_2_name = basename(fastq_file_2, "_2.fastq") + "_downsample_~{final_depth}x"

  output {
    File downsample_file_1 = fastq_downsample_1_name + ".1.fastq"
    File downsample_file_2 = fastq_downsample_2_name + ".2.fastq"
  }

  command <<<
    set -euo pipefail

    #counts the number of reads in each fastq file and prints for qc
    count1=$(bash /opt/count_fastq.sh ~{fastq_file_1})
    echo ${count1} 

    count2=$(bash /opt/count_fastq.sh ~{fastq_file_2})
    echo ${count2}

    #calculates half of the initial and final depth for the number of reads to sample in each fq file
    initial=$(echo ~{start_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
    final=$(echo ~{final_depth} | awk ' { printf "%0.2f\n", ($1 / 2); } ')
    echo ${initial}
    echo ${final}

    if [ $count1 -eq $count2 ] #checks that the paired-end fastq files have the same number of reads
      then
      #half of the total read coverage (i.e. 15x for 30x)
        freads=$(bash /opt/calcRD.sh ${initial} ${final} ${count1})  

        echo "freads: ${freads}"

        #uses now the calculated freads (1/2 of final read coverage wanted) in each fastq file
        fastq-sample -n ${freads} --seed ~{seed} -o ~{fastq_downsample_1_name} ~{fastq_file_1} ~{fastq_file_2}
        #output name will be: sample_id.final_downsample_~{final_depth}x.[1/2].fastq

      else
        echo "Error: counts don't match up!"
        echo ${count1}
        echo ${count2}
    fi
  >>>

  runtime {
    docker: downsample_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}


