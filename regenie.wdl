version 1.0

## Version 08-19-2021
##
## This WDL workflow runs Regenie.
## This workflow assumes users have thoroughly read the Regenie docs for caveats and details.
## Regenie's documentation: https://rgcgithub.github.io/regenie/options/
##
##  Step1 - Made with BGEN files in mind to use as the input genetic data file. - (BGEN version 1.2, 8-bit probabilities).
##  Step2 - Use separate .bed, .bim, .fam files for each chromosome. (If testing chr 1-22, there should be 22 separate files).
##  PLINK can be used to convert bed, bim, and fam files to BGEN files outside of this workflow.
##  PLINK can also be used to convert pgen, pvar, and psam files to BGEN files outside of this workflow.
##
## Cromwell version support - Successfully tested on v66
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow Regenie {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        File bgen_step1
        File bed_files_step2_file # File containing tab-separated list of: .bed file, .bim file, .fam file, and corresponding chromosome number - for all files used for Step2.
        String fit_bin_out_name = "fit_bin_out" # File prefix for the list of predictions produced in Step1.
        File? fit_bin_out
        File? sample
        File? keep
        File? remove
        File? extract
        File? exclude
        File? phenoFile
        String? phenoCol
        String? phenoColList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList
        File? pred
        Array[Int] chr_list # List of chromosomes used for analysis.
        Array[String] phenotype_names # Phenotypes you want to analyze. (Column names).
        String regenie_docker = "briansha/regenie:v2.2.2_boost" # Compiled with Boost IOSTREAM: https://github.com/rgcgithub/regenie/wiki/Using-docker
        String r_base_docker = "briansha/regenie_r_base:4.1.0"  # Ubuntu 18.04, R 4.1.0, and a few Ubuntu and R packages.
    }
    Array[Array[String]] bed_files_step2 = read_tsv(bed_files_step2_file)

    call RegenieStep1WholeGenomeModel {
        input:
            bgen_step1 = bgen_step1,
            keep = keep,
            fit_bin_out_name = fit_bin_out_name,
            covarFile = covarFile,
            exclude = exclude,
            phenoFile = phenoFile,
            remove = remove,
            docker = regenie_docker
    }

    scatter (files in bed_files_step2) {
      call RegenieStep2AssociationTesting {
          input:
              bed_step2 = files[0],
              bim_step2 = files[1],
              fam_step2 = files[2],
              chr_name = files[3],
              keep = keep,
              covarFile = covarFile,
              exclude = exclude,
              phenoFile = phenoFile,
              remove = remove,
              pred = RegenieStep1WholeGenomeModel.fit_bin_out,
              docker = regenie_docker,
              output_locos = RegenieStep1WholeGenomeModel.output_locos
      }
    }

    call join_Output {
      input:
        output_files = RegenieStep2AssociationTesting.test_bin_out_firth, # Refers implicitly to the entire array of files that were scattered.
        chr_list = chr_list,
        phenotype_names = phenotype_names,
        docker = r_base_docker
    }

    call Plots {
      input:
        chr_list = chr_list,
        docker = r_base_docker,
        file_input = join_Output.outputs
    }

    output {
         Array[File] output_plots = Plots.output_plots
         Array[File] output_regenie = Plots.output_regenie
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow runs Regenie - see the README on the github for more information - https://github.com/briansha/Regenie_WDL."
    }
}

# In Step 1, the whole genome regression model is fit to the traits
# and a set of genomic predictions are produced as output.
task RegenieStep1WholeGenomeModel {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        String docker
        String fit_bin_out_name # File prefix for the list of predictions produced.
        File bgen_step1
        Int bsize = 1000  # 1000 is recommended by regenie's documentation
        Float memory = 3.5
        Int? disk_size_override
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0

        File? sample
        File? keep
        File? remove
        File? extract
        File? exclude
        File? phenoFile
        String? phenoCol
        String? phenoColList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList
        File? pred

        Boolean ref_first = false
        Boolean lowmem = false
        Boolean lowmem_prefix_flag = false
        String? lowmem_prefix = false
        Boolean bt = false # Specifies both phenotype file and covariate files contain binary values.
        Boolean cc12 = false
        Int? cv
        Boolean loocv = false
        String? split_l0
        File? run_l0
        File? run_l1
        Boolean keep_l0 = false
        Boolean print_prs = false
        Boolean force_step1 = false
        Int? nb
        Boolean strict = false
        Boolean ignore_pred = false
        Boolean use_prs = false
        Boolean gz = false
        Boolean force_impute = false
        Boolean write_samples = false
        Boolean print_pheno = false
        Boolean firth = false
        Boolean approx = false
        Boolean firth_se = false
        Boolean spa = false
        Float? pThresh
        String? test
        Int? chr
        String? chrList
        String? range
        Float? minMAC
        Float? minINFO
        Int? nauto
        Int? maxCatLevels
        Int? niter
        Int? maxstep_null
        Int? maxiter_null
        Boolean debug = false
        Boolean verbose = false
        Boolean help = false
    }
    Float genotype_size = size(bgen_step1, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * genotype_size)])

    command <<<
        set -euo pipefail
        regenie \
        --step 1 \
        --bgen=~{bgen_step1} \
        --phenoFile=~{phenoFile} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(sample) then "--sample=~{sample} " else " "} \
        ~{if ref_first then "--ref-first " else " "} \
        ~{if defined(keep) then "--keep=~{keep} " else " "} \
        ~{if defined(remove) then "--remove=~{remove} " else " "} \
        ~{if defined(extract) then "--extract=~{extract} " else " "} \
        ~{if defined(exclude) then "--exclude=~{exclude} " else " "} \
        ~{if defined(phenoCol) then "--phenoCol=~{phenoCol} " else " "} \
        ~{if defined(phenoColList) then "--phenoColList=~{phenoColList} " else " "} \
        ~{if defined(covarCol) then "--covarCol=~{covarCol} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        ~{if defined(pred) then "--pred=~{pred} " else " "} \
        --bsize=~{bsize} \
        ~{if lowmem then "--lowmem " else " "} \
        ~{if lowmem_prefix_flag then "--lowmem-prefix=~{lowmem_prefix} " else " "} \
        ~{if bt then "--bt " else " "} \
        ~{if cc12 then "-1,--c12 " else " "} \
        ~{if defined(cv) then "--cv=~{cv} " else " "} \
        ~{if loocv then "--loocv " else " "} \
        ~{if defined(split_l0) then "--split_l0=~{split_l0} " else " "} \
        ~{if defined(run_l0) then "--run_l0=~{run_l0} " else " "} \
        ~{if defined(run_l1) then "--run_l1=~{run_l1} " else " "} \
        ~{if keep_l0 then "--keep_l0 " else " "} \
        ~{if print_prs then "--print_prs " else " "} \
        ~{if force_step1 then "--force_step1 " else " "} \
        ~{if defined(nb) then "--nb=~{nb} " else " "} \
        ~{if strict then "--strict " else " "} \
        ~{if ignore_pred then "--ignore_pred " else " "} \
        ~{if use_prs then "--use_prs " else " "} \
        ~{if gz then "--gz " else " "} \
        ~{if force_impute then "--force_impute " else " "} \
        ~{if write_samples then "--write_samples " else " "} \
        ~{if print_pheno then "--print_pheno " else " "} \
        ~{if firth then "--firth " else " "} \
        ~{if approx then "--approx " else " "} \
        ~{if firth_se then "--firth_se " else " "} \
        ~{if spa then "--spa " else " "} \
        ~{if defined(pThresh) then "--pThresh=~{pThresh} " else " "} \
        ~{if defined(test) then "--test=~{test} " else " "} \
        ~{if defined(chr) then "--chr=~{chr} " else " "} \
        ~{if defined(chrList) then "--chrList=~{chrList} " else " "} \
        ~{if defined(range) then "--range=~{range} " else " "} \
        ~{if defined(minMAC) then "--minMAC=~{minMAC} " else " "} \
        ~{if defined(minINFO) then "--minINFO=~{minINFO} " else " "} \
        ~{if defined(nauto) then "--nauto=~{nauto} " else " "} \
        ~{if defined(maxCatLevels) then "--maxCatLevels=~{maxCatLevels} " else " "} \
        ~{if defined(niter) then "--niter=~{niter} " else " "} \
        ~{if defined(maxstep_null) then "--maxstep_null=~{maxstep_null} " else " "} \
        ~{if defined(maxiter_null) then "--maxiter_null=~{maxiter_null} " else " "} \
        ~{if debug then "--debug " else " "} \
        ~{if verbose then "--verbose " else " "} \
        ~{if help then "--help " else " "} \
        --out fit_bin_out
    >>>

    output {
        File fit_bin_out = "${fit_bin_out_name}_pred.list" # Refers to the list of loco files written. Lists the files as being located in the current working directory.
        Array[File] output_locos = glob("*.loco") # Writes n loco files for n phenotypes.
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# In Step 2, a set of imputed SNPs are tested for association using a Firth logistic regression model.
task RegenieStep2AssociationTesting {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input {
        String docker
        File bed_step2
        File bim_step2
        File fam_step2
        Array[File] output_locos
        File? fit_bin_out
        Float memory = 3.5
        Int? disk_size_override
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0

        Int bsize = 1000  # 1000 is recommended by regenie's documentation
        File? pred # File containing predictions from Step 1
        String? chr_name
        File? sample
        File? keep
        File? remove
        File? extract
        File? exclude
        File? phenoFile
        String? phenoCol
        String? phenoColList
        File? covarFile
        String? covarCol
        String? covarColList
        String? catCovarList

        Boolean ref_first = false
        Boolean lowmem = false
        Boolean lowmem_prefix_flag = false
        String? lowmem_prefix
        Boolean bt = false
        Boolean cc12 = false
        Int? cv
        Boolean loocv = false
        String? split_l0
        File? run_l0
        File? run_l1
        Boolean keep_l0 = false
        Boolean print_prs = false
        Int? nb
        Boolean strict = false
        Boolean ignore_pred = false
        Boolean use_prs = false
        Boolean gz = false
        Boolean force_impute = false
        Boolean write_samples = false
        Boolean print_pheno = false
        Boolean firth = false
        Boolean approx = false
        Boolean firth_se = false
        Boolean spa = false
        Float? pThresh
        String? test
        Int? chr
        String? chrList
        String? range
        Float? minMAC
        Float? minINFO
        Int? nauto
        Int? maxCatLevels
        Int? niter
        Int? maxstep_null
        Int? maxiter_null
        Boolean debug = false
        Boolean verbose = false
        Boolean help = false

        File? anno_file
        File? set_list
        File? extract_sets
        File? exclude_sets
        String? extract_setlist
        String? exclude_setlist
        File? aaf_file
        File? mask_def
        Float? aaf_bins
        String? build_mask
        Boolean singleton_carrier = false
        Boolean write_mask = false
        Boolean skip_test = false
        String? mask_lovo
        Boolean check_burden_files = false
        Boolean strict_check_burden = false
    }
    Float dosage_size = size(bed_step2, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * dosage_size)])

    # Loco files are moved to the current working directory due to the list of predictions (pred) expecting them to be there.
    command <<<
        set -euo pipefail
        for file in ~{sep=' ' output_locos}; do \
          mv $file .; \
        done
        regenie \
        --step 2 \
        --bed=~{sub(bed_step2,'\\.bed$','')} \
        --phenoFile=~{phenoFile} \
        ~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} \
        ~{if defined(sample) then "--sample=~{sample} " else " "} \
        ~{if ref_first then "--ref-first " else " "} \
        ~{if defined(keep) then "--keep=~{keep} " else " "} \
        ~{if defined(remove) then "--remove=~{remove} " else " "} \
        ~{if defined(extract) then "--extract=~{extract} " else " "} \
        ~{if defined(exclude) then "--exclude=~{exclude} " else " "} \
        ~{if defined(phenoCol) then "--phenoCol=~{phenoCol} " else " "} \
        ~{if defined(phenoColList) then "--phenoColList=~{phenoColList} " else " "} \
        ~{if defined(covarCol) then "--covarCol=~{covarCol} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(catCovarList) then "--catCovarList=~{catCovarList} " else " "} \
        ~{if defined(pred) then "--pred=~{pred} " else " "} \
        --bsize=~{bsize} \
        ~{if lowmem then "--lowmem " else " "} \
        ~{if lowmem_prefix_flag then "--lowmem-prefix=~{lowmem_prefix} " else " "} \
        ~{if bt then "--bt " else " "} \
        ~{if cc12 then "-1,--c12 " else " "} \
        ~{if defined(cv) then "--cv=~{cv} " else " "} \
        ~{if loocv then "--loocv " else " "} \
        ~{if defined(split_l0) then "--split_l0=~{split_l0} " else " "} \
        ~{if defined(run_l0) then "--run_l0=~{run_l0} " else " "} \
        ~{if defined(run_l1) then "--run_l1=~{run_l1} " else " "} \
        ~{if keep_l0 then "--keep_l0 " else " "} \
        ~{if print_prs then "--print_prs " else " "} \
        ~{if defined(nb) then "--nb=~{nb} " else " "} \
        ~{if strict then "--strict " else " "} \
        ~{if ignore_pred then "--ignore_pred " else " "} \
        ~{if use_prs then "--use_prs " else " "} \
        ~{if gz then "--gz " else " "} \
        ~{if force_impute then "--force_impute " else " "} \
        ~{if write_samples then "--write_samples " else " "} \
        ~{if print_pheno then "--print_pheno " else " "} \
        ~{if firth then "--firth " else " "} \
        ~{if approx then "--approx " else " "} \
        ~{if firth_se then "--firth_se " else " "} \
        ~{if spa then "--spa " else " "} \
        ~{if defined(pThresh) then "--pThresh=~{pThresh} " else " "} \
        ~{if defined(test) then "--test=~{test} " else " "} \
        ~{if defined(chr) then "--chr=~{chr} " else " "} \
        ~{if defined(chrList) then "--chrList=~{chrList} " else " "} \
        ~{if defined(range) then "--range=~{range} " else " "} \
        ~{if defined(minMAC) then "--minMAC=~{minMAC} " else " "} \
        ~{if defined(minINFO) then "--minINFO=~{minINFO} " else " "} \
        ~{if defined(nauto) then "--nauto=~{nauto} " else " "} \
        ~{if defined(maxCatLevels) then "--maxCatLevels=~{maxCatLevels} " else " "} \
        ~{if defined(niter) then "--niter=~{niter} " else " "} \
        ~{if defined(maxstep_null) then "--maxstep_null=~{maxstep_null} " else " "} \
        ~{if defined(maxiter_null) then "--maxiter_null=~{maxiter_null} " else " "} \
        ~{if debug then "--debug " else " "} \
        ~{if verbose then "--verbose " else " "} \
        ~{if help then "--help " else " "} \
        ~{if defined(anno_file) then "--anno_file=~{anno_file} " else " "} \
        ~{if defined(set_list) then "--set_list=~{set_list} " else " "} \
        ~{if defined(extract_sets) then "--extract_sets=~{extract_sets} " else " "} \
        ~{if defined(exclude_sets) then "--exclude_sets=~{exclude_sets} " else " "} \
        ~{if defined(extract_setlist) then "--extract_setlist=~{extract_setlist} " else " "} \
        ~{if defined(exclude_setlist) then "--exclude_setlist=~{exclude_setlist} " else " "} \
        ~{if defined(aaf_file) then "--aaf_file=~{aaf_file} " else " "} \
        ~{if defined(mask_def) then "--mask_def=~{mask_def} " else " "} \
        ~{if defined(aaf_bins) then "--aaf_bins=~{aaf_bins} " else " "} \
        ~{if defined(build_mask) then "--build_mask=~{build_mask} " else " "} \
        ~{if singleton_carrier then "--singleton_carrier " else " "} \
        ~{if write_mask then "--write_mask " else " "} \
        ~{if skip_test then "--skip_test " else " "} \
        ~{if defined(mask_lovo) then "--mask_lovo=~{mask_lovo} " else " "} \
        ~{if check_burden_files then "--check_burden_files " else " "} \
        ~{if strict_check_burden then "--strict_check_burden " else " "} \
        --out ~{chr_name}
    >>>

    output {
        Array[File] test_bin_out_firth = glob("*.regenie")
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Join together all output relating to the scattered tasks in Step 2.
# For each phenotype, one file will be made containing all analysis from Step 1 on each chromosome used in testing.
task join_Output {
  input {
    Array[Array[File]] output_files # All output from Step 2.
    Array[String] phenotype_names
    Array[Int] chr_list
    String docker
    Float memory = 3.5
    Int? disk_size_override
    Int cpu = 1
    Int preemptible = 1
    Int maxRetries = 0
  }
  Array[File] all_output_files = flatten(output_files)
  Float regenie_files_size = size(all_output_files, "GiB")
  Int disk = select_first([disk_size_override, ceil(30.0 + 2.0 * regenie_files_size)])

  command <<<
        set -euo pipefail
        for array in ~{sep=' ' all_output_files}; do \
          for file in $array; do \
            mv $file .; \
        done done
        for pheno in ~{sep=' ' phenotype_names}; do \
          for chr in ~{sep= ' ' chr_list}; do \
            cat ${chr}_${pheno}.regenie >> $pheno.regenie; \
            rm ${chr}_${pheno}.regenie; \
        done done
  >>>

  output {
    Array[File] outputs = glob("*.regenie")
  }

  runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
  }
}

# QQ and Manhattan plots
task Plots {
  input {
    Array[String] phenotype_names_regenie # Format: "<phenotype_name>.regenie". Names of the files produced in join_Output.
    Array[Int] chr_list
    Array[File] file_input
    String docker
    Float memory = 16.0
    Int? disk_size_override
    Int cpu = 1
    Int preemptible = 1
    Int maxRetries = 0
  }
  Float regenie_files_size = size(file_input, "GiB")
  Int disk = select_first([disk_size_override, ceil(10.0 + regenie_files_size)])

  # Plots are produced for each phenotype.
  # For each phenotype, a file containing all of the hits from Step 2 is output.
  # For each phenotype, a file containing a subset of all of the hits where "-LOG10P > 1.3" from Step 2 is output.
  command <<<
    set -euo pipefail
    for file in ~{sep=' ' file_input}; do \
      awk '$12 > 1.3' $file >> ${file%.regenie}_subset.regenie; \
      mv ${file%.regenie}_subset.regenie .; \
      mv $file .; \
    done
    R --no-save --args ~{sep=' ' phenotype_names_regenie} <<RSCRIPT
    library(data.table)
    library(qqman)
    args <- commandArgs(trailingOnly = TRUE)
    for (file in args) {
      regenie_output <- fread(file)
        regenie_ADD_subset <-subset.data.frame(regenie_output, TEST=="ADD")
        regenie_ADD_subset[,"CHROM"] <-as.numeric(unlist(regenie_ADD_subset[,"CHROM"]))
        regenie_ADD_subset[,"LOG10P"] <-as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))
        regenie_ADD_subset[,"GENPOS"] <-as.numeric(unlist(regenie_ADD_subset[,"GENPOS"]))
        qq_plot = substr(file,1,nchar(file)-8)
        qq_plot = paste(qq_plot, "qqplot.png", sep="_")
        png(qq_plot)
        p = 10 ^ (-1 * (as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))))
        print(qq(p))
        dev.off()
        manhattan_plot = substr(file,1,nchar(file)-8)
        manhattan_plot = paste(manhattan_plot, "manhattan.png", sep="_")
        png(manhattan_plot)
        print(manhattan(regenie_ADD_subset, chr="CHROM", bp="GENPOS", snp="ID", p="LOG10P", logp=FALSE, annotatePval = 1E-5))
        dev.off()
      }
    RSCRIPT
  >>>

  output {
    Array[File] output_plots = glob("*.png")
    Array[File] output_regenie = glob("*.regenie")
  }

  runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
  }
}
