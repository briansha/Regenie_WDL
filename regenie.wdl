version development

## Version 04-05-2021
##
## This WDL workflow runs regenie.
## Made with BGEN files in mind to use as the input genetic data file. - (BGEN version 1.2, 8-bit probabilities)
## PLINK can be used to convert bed, bim, and fam files to BGEN files outside of this workflow.
## PLINK can also be used to convert pgen, pvar, and psam files to BGEN files outside of this workflow.
##
## Cromwell version support - Successfully tested on v59
##
## Distributed under terms of the MIT License
## Copyright (c) 2020-2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow Regenie {

    input {
        File bgen_step1
        File bgen_step2
        String fit_bin_out_name
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
        Int? memory
        Int? disk
        Int? threads
        String docker_image
        String docker_image_R = "r-base:4.0.3"
        Array[Int] chr_list # List of chromosomes for Step 2.
        Int? n_phenotypes # Number of phenotypes used for Step 1.
        Array[String] phenotype_names
        
    }

  
    call Step1 { 
        input: 
            bgen_step1 = bgen_step1,
            keep = keep,
            fit_bin_out_name = fit_bin_out_name,
            covarFile = covarFile,
            exclude = exclude,
            phenoFile = phenoFile,
            remove = remove,
            memory = memory,
            disk = disk,
            threads = threads,
            docker_image = docker_image
    }
    
    scatter (chromosome in chr_list) {
      call Step2 { 
          input: 
              bgen_step2 = bgen_step2,
              keep = keep,
              chr = chromosome, # May need to change to chrList = chr_list later, as chr is Int only...not accounting for X, Y, etc.
              covarFile = covarFile,
              exclude = exclude,
              phenoFile = phenoFile,
              remove = remove,
              pred = Step1.fit_bin_out,
              memory = memory,
              disk = disk,
              threads = threads,
              docker_image = docker_image,
              output_locos = Step1.output_locos
      }
    }
    
    call join_Output {
      input:
        output_files = Step2.test_bin_out_firth, # Refers implicitly to the entire array of files that were scattered.
        chr_list = chr_list,
        memory = memory,
        disk = disk,
        threads = threads,
        phenotype_names = phenotype_names,
        docker_image_R = docker_image_R
    }
    
    call Plots {
      input:
        chr_list = chr_list,
        memory = memory,
        disk = disk,
        threads = threads,
        docker_image_R = docker_image_R,
        file_input = join_Output.outputs
    }

    output {
         Array[File] output_plots = Plots.output_plots
         Array[File] output_regenie = Plots.output_regenie
    }
    
    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "Run REGENIE"
       }
}

task Step1 {

    input {
        String docker_image
        String fit_bin_out_name
        File bgen_step1
        Int bsize  # 1000 is recommended by regenie's documentation
        Int? memory
        Int? disk
        Int? threads
        
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
        
        Boolean ref_first
        Boolean lowmem
        Boolean lowmem_prefix_flag
        String? lowmem_prefix
        Boolean bt  # Specifies both phenotype file and covariate files contain binary values.
        Boolean cc12
        Int? cv
        Boolean loocv
        String? split_l0
        File? run_l0
        File? run_l1
        Boolean keep_l0
        Boolean print_prs
        Boolean force_step1
        Int? nb
        Boolean strict
        Boolean ignore_pred
        Boolean use_prs
        Boolean gz
        Boolean force_impute
        Boolean write_samples
        Boolean print_pheno
        Boolean firth
        Boolean approx
        Boolean firth_se
        Boolean spa
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
        Boolean debug
        Boolean verbose
        Boolean help
    }

    command <<<
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
        File fit_bin_out = "${fit_bin_out_name}_pred.list"
        Array[File] output_locos = glob("*.loco") # Writes n_files for n_phenotypes.
    }

    runtime {
        docker: docker_image
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: threads
    }
}

task Step2 {

    input {
        String docker_image
        File bgen_step2
        Array[File] output_locos
        File? fit_bin_out
        Int? memory
        Int? disk
        Int? threads
        Int bsize  # 1000 is recommended by regenie's documentation
        
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
        
        Boolean ref_first
        Boolean lowmem
        Boolean lowmem_prefix_flag
        String? lowmem_prefix
        Boolean bt
        Boolean cc12
        Int? cv
        Boolean loocv
        String? split_l0
        File? run_l0
        File? run_l1
        Boolean keep_l0
        Boolean print_prs
        Int? nb
        Boolean strict
        Boolean ignore_pred
        Boolean use_prs
        Boolean gz
        Boolean force_impute
        Boolean write_samples
        Boolean print_pheno
        Boolean firth
        Boolean approx
        Boolean firth_se
        Boolean spa
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
        Boolean debug
        Boolean verbose
        Boolean help
        
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
        Boolean singleton_carrier
        Boolean write_mask
        Boolean skip_test
        String? mask_lovo
        Boolean check_burden_files
        Boolean strict_check_burden
    }


    command <<<
        for file in ~{sep=' ' output_locos}; do \
          mv $file .; \
        done
        regenie \
        --step 2 \
        --bgen=~{bgen_step2} \
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
        --out ~{chr}
    >>>
    
    output {
        Array[File] test_bin_out_firth = glob("*.regenie")
    }

    runtime {
        docker: docker_image
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: threads
    }
}

task join_Output {

  input {
    Array[Array[File]] output_files
    Array[String] phenotype_names
    Array[Int] chr_list
    Int? memory
    Int? disk
    Int? threads
    String docker_image_R
  }
  
  Array[File] all_output_files = flatten(output_files)
  
  command <<<
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
        docker: docker_image_R
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: threads
  }
}
  task Plots {

  input {
    Array[String] phenotype_names_regenie
    Array[Int] chr_list
    Int? memory
    Int? disk
    Int? threads
    String docker_image_R
    Array[File] file_input
  }
  
  command <<<
    for file in ~{sep=' ' file_input}; do \
      awk '$13 > 1.3' $file >> ${file%.regenie}_subset.regenie; \
      mv ${file%.regenie}_subset.regenie .; \
      mv $file .; \
    done
    R --no-save --args ~{sep=' ' phenotype_names_regenie} <<RSCRIPT
    install.packages("data.table", dependencies=TRUE)
    library(data.table)
    install.packages("qqman", dependencies=TRUE)
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
      print(qq(as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))))
      dev.off()
      manhattan_plot = substr(file,1,nchar(file)-8)
      manhattan_plot = paste(manhattan_plot, "manhattan.png", sep="_")
      png(manhattan_plot)
      print(manhattan(regenie_ADD_subset, chr="CHROM", bp="GENPOS", snp="ID", p="LOG10P", annotatePval = 1E-5))
      dev.off()
      }
    RSCRIPT
  >>>
  
  output {
    Array[File] output_plots = glob("*.png")
    Array[File] output_regenie = glob("*.regenie")
  }
  
  runtime {
        docker: docker_image_R
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: threads
  }
}
