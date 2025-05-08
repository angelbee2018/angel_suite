# NCI PARALLEL FUNCTIONS
## skibidi toilet if you run this script your computer will explode lmfao
## R wrapper for launching child processes using the mpirun ... nci-parallel system in NCI Gadi

# nciparallel_pmap
## it allows purrr functions to be done in parallel across nodes
## how it works: this is the final form. it is ALL disk-based.
## 1. triaging and chunking
## 2. autoscoping and export to disk
###   - save .l in chunks
###   - inject prep code into .f and write to file
###   - export environment variables as usual
## 3. worker commissioning
###   - Launch: system command
###   - Process poll: ????
###   - Status updates: ????
###   - Completion poll: read status from disk
###   - Process killing: ????
## 4. termination
###   - splicing: is always ordered
###   - save as chunks is NOT needed because if the result was so gigantic that it couldnt be spliced into a single instance, the user code should already account for that.

# IMPORTANT ###
## It is essential that you run the following commands in bash BEFORE launching your compute job. 
## This adds your ssh key to yourself so that the polling algorithm can fetch the process status without being stalled by password requests.
# # If you don’t already have one—create a key with no passphrase:
# ssh-keygen -t ed25519 -N '' -f ~/.ssh/id_ed25519
# cat ~/.ssh/id_ed25519.pub >>  ~/.ssh/authorized_keys
# chmod u=rw,go= ~/.ssh/authorized_keys
###############

# KNOWN BUGS THAT ARE OUT OF OUR CONTROL ###
# - if running in ARE, you need to "File > Quit Session" before each run, otherwise `echo $PBS_JOBID` will reflect an old job
############################################

## .dedicated_mgmt_rank: nci-parallel has a feature: '--dedicated,-d dedicate rank 0 to task management'
## .run_in_background: frees up the R process to do other things. but by definition, it will NEVER return a result directly to R and chunks are always saved to disk.
nciparallel_pmap <- function(
    .l, .f, 
    .no_workers = 1, .no_chunks = 1, .no_workers_per_node = 1, .no_cores_per_worker = 1, .r_process_spawn_wait_time = 120, .commission_mode = c("nci_parallel_managed_round_robin", "r_managed_round_robin"), .dedicated_mgmt_rank = FALSE, .use_parallel_for_mgmt = FALSE, .run_in_background = FALSE,
    .job_name = NULL, .patience = 5,
    .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL, 
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE, 
    .status_messages_dir = NULL, .progress = TRUE,
    .debug = FALSE) {
  
  # DEBUG ###
  # .l = list(
  #   "a1" = list.files(path = vcf_dir, pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$") %>% gsub(pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$", replacement = ""),
  #   "a2" = 1:length(list.files(path = vcf_dir, pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$"))
  # )
  # .no_workers = 16
  # .no_chunks = length(list.files(path = vcf_dir, pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$"))
  # .no_workers_per_node = 8
  # .no_cores_per_worker = 28
  # .r_process_spawn_wait_time = 120
  # .commission_mode = "nci_parallel_managed_round_robin"
  # .f = function(a1, a2) {
  # 
  #   # DEBUG ###
  #   # a1 <- list.files(path = vcf_dir, pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$") %>% gsub(pattern = "\\_OM\\_DNA\\_blood\\.hard\\-filtered\\.vcf$", replacement = "") %>% .[[4]]
  #   ###########
  # 
  #   cat(paste("Progress: ", a2, "\n", sep = ""))
  # 
  #   vcfr_blood <- vcfR::read.vcfR(file = paste(vcf_dir, a1, "_OM_DNA_blood.hard-filtered.vcf", sep = ""))
  # 
  #   tibble_fixed_data_reduced <- vcfr_blood@fix %>% tibble::as_tibble()
  #   tibble_fixed_data_reduced <- tibble_fixed_data_reduced[tibble_fixed_data_reduced$FILTER == "PASS", ]
  #   tibble_fixed_data_reduced <- tibble_fixed_data_reduced[, c("CHROM", "POS", "REF", "ALT", "INFO")]
  #   colnames(tibble_fixed_data_reduced)[colnames(tibble_fixed_data_reduced) == "INFO"] <- paste("AC_", a1, sep = "")
  #   tibble_fixed_data_reduced[grep(x = tibble_fixed_data_reduced[[paste("AC_", a1, sep = "")]], pattern = "AC\\=", invert = TRUE), paste("AC_", a1, sep = "")] <- "1"
  #   tibble_fixed_data_reduced[, paste("AC_", a1, sep = "")] <- tibble_fixed_data_reduced[[paste("AC_", a1, sep = "")]] %>% gsub(pattern = "AC\\=([^;]+);.*", replacement = "\\1")
  #   tibble_fixed_data_reduced <- split_delimited_columns_in_table2(input_table = tibble_fixed_data_reduced, target_colnames = c("ALT", paste("AC_", a1, sep = "")), split = "\\,")
  # 
  #   # tibble_genotype_data <- vcfr_blood@gt %>% tibble::as_tibble()
  # 
  #   L1_tibble_matrix_recipient_variant_presence <- dplyr::full_join(tibble_variants_database, tibble_fixed_data_reduced)
  #   L1_tibble_matrix_recipient_variant_presence[is.na(L1_tibble_matrix_recipient_variant_presence[[paste("AC_", a1, sep = "")]]), paste("AC_", a1, sep = "")] <- 0
  # 
  #   # data.table::fwrite(x = L1_tibble_matrix_recipient_variant_presence[, paste("AC_", a1, sep = "")], file = paste(r_results_dir, paste(vector_experiment_tag, collapse = "_"), "_vector_variant_presence_recipient_temp_", a2, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  # 
  #   return(L1_tibble_matrix_recipient_variant_presence[[paste("AC_", a1, sep = "")]])
  # 
  #   # return(NULL)
  # 
  # }
  # .dedicated_mgmt_rank = FALSE
  # .use_parallel_for_mgmt = FALSE
  # .run_in_background = FALSE
  # .job_name = NULL
  # .patience = 5
  # .globals_save_compress = TRUE
  # .re_export = TRUE
  # .globals_mode = "auto"
  # .user_global_objects = NULL
  # .intermediate_files_dir = NULL
  # .keep_intermediate_files = FALSE
  # .status_messages_dir = NULL
  # .progress = TRUE
  # .debug = FALSE
  ###########
  
  for (i in c("globals", "parallel", "utils", "lubridate", "qs")) {
    if (require(i, character.only = TRUE) == FALSE) {
      stop(paste("Package \"", i, "\" not found. Please install using `install.packages` or `BiocManager::install`", sep = ""))
    }
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    if (type.convert(system(command = "ulimit -n", intern = TRUE), as.is = TRUE) < 65536) {
      warning(paste("System max. open file limit is less than the recommended 65536. Current limit is set to: ", system(command = "ulimit -n", intern = TRUE), ". To fix this, please re-run R from bash terminal after having set `ulimit -n 65536` to avoid possible errors with large jobs and/or large number of workers/chunks.", sep = ""))
    }
  }
  
  .parent_pid <- Sys.getpid()
  
  .epoch_time <- as.numeric(Sys.time())*1E9
  .temp_dir <- paste(tempdir(), "_", .epoch_time, "/", sep = "")
  
  if (is.null(.job_name)) {
    .job_name <- paste("nciparallel_pmap_", .epoch_time, sep = "")
  }
  
  message(paste("Job name: ", .job_name, sep = ""))
  
  if (is.null(.intermediate_files_dir)) {
    .intermediate_files_dir <- .temp_dir
  }
  if (!dir.exists(.intermediate_files_dir)) {
    dir.create(.intermediate_files_dir, recursive = TRUE)
  }
  
  if (.debug == TRUE) {
    message(paste(".keep_intermediate_files: ", .keep_intermediate_files, sep = ""))
  }
  
  message(paste("Intermediate files will be saved to disk at: ", .intermediate_files_dir, sep = ""))
  
  if (is.null(.status_messages_dir)) {
    .status_messages_dir <- .temp_dir
  }
  if (!dir.exists(.status_messages_dir)) {
    dir.create(.status_messages_dir, recursive = TRUE)
  }
  
  message(paste("Status messages will be saved to disk at: ", .status_messages_dir, sep = ""))
  
  # SCOPING ####
  
  vector_global_packages <- globals::packagesOf(globals::globalsOf(.f, mustExist = FALSE)) %>% setdiff(., c("base", "rlang"))
  
  # write vector_global_packages
  qs::qsave(x = vector_global_packages, file = paste(.intermediate_files_dir, "/", .job_name, "_vector_global_packages.qs", sep = ""))
  
  if (.debug == TRUE) {
    print("vector_global_packages at the surface level of the function")
    print(vector_global_packages)
  }
  
  # save the current env into temp path
  if (.re_export == FALSE & dir.exists(.intermediate_files_dir)) {
    
    message("Temp object data file found. Will load that instead of re-exporting")
    
  } else {
    
    if (.globals_mode == "global") {
      message("Writing globals to disk")
      save.image(file = .intermediate_files_dir, compress = .globals_save_compress)
    } else if (.globals_mode == "auto") {
      
      vector_global_variables <- c(globals::findGlobals(.f), .user_global_objects)
      
      if (.debug == TRUE) {
        print(vector_global_variables)
      }
      
      temp_global_variables_uncollected <- setdiff(vector_global_variables, ls(all.names = TRUE))
      
      # recursively traverse the frame stack from the bottom up until we collect everything
      temp_current_frame <- sys.nframe()
      
      while (length(temp_global_variables_uncollected) > 0 & temp_current_frame > -1) {
        
        for (j in temp_global_variables_uncollected) {
          
          assign(x = j, value = dynGet(x = j, ifnotfound = NULL, minframe = temp_current_frame))
          
          if (.debug == TRUE) {
            message("temp_current_frame: ", temp_current_frame)
            message(j)
            print(get(x = j))
          }
          
          if (is.null(get(x = j)) == TRUE) {
            rm(list = j)
          }
          
        }
        
        temp_global_variables_uncollected <- setdiff(temp_global_variables_uncollected, ls(all.names = TRUE))
        
        temp_current_frame <- temp_current_frame - 1
        
      }
      
      # deal with the special case of nested function definitions
      vector_all_function_names_in_environment <- lsf.str() %>% as.character
      
      temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- vector_all_function_names_in_environment
      temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- character()
      
      # recursively inspect functions until there are no more functions within the functions
      # nested packages -> added to global vector of package names
      # nested functions -> copied into environment using `dynGet` + uncollected function names are added to the uncollected list + added to vector_global_variables + flagged for further internal inspection
      while (length(temp_vector_names_of_newly_discovered_nested_functions_for_inspection) > 0) {
        
        for (temp_function_to_inspect in purrr::map(.x = temp_vector_names_of_newly_discovered_nested_functions_for_inspection, .f = ~get(.x))) {
          
          vector_global_packages <- c(vector_global_packages, globals::packagesOf(globals::globalsOf(temp_function_to_inspect, mustExist = FALSE)) %>% setdiff(., c("base", "rlang"))) %>% unique
          
          temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- c(temp_vector_names_of_newly_discovered_nested_functions_for_inspection0, setdiff(globals::findGlobals(temp_function_to_inspect), temp_global_variables_uncollected))
          
        }
        
        # flag for further internal inspection
        temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- temp_vector_names_of_newly_discovered_nested_functions_for_inspection0
        
        # copy into env using dynGet
        # recursively traverse the frame stack from the bottom up until we collect everything
        temp_current_frame <- sys.nframe()
        
        while (length(temp_vector_names_of_newly_discovered_nested_functions_for_inspection0) > 0 & temp_current_frame > -1) {
          
          for (j in temp_vector_names_of_newly_discovered_nested_functions_for_inspection0) {
            
            assign(x = j, value = dynGet(x = j, ifnotfound = NULL, minframe = temp_current_frame))
            
            if (.debug == TRUE) {
              message("temp_current_frame: ", temp_current_frame)
              message(j)
              print(get(x = j))
            }
            
            if (is.null(get(x = j)) == TRUE) {
              rm(list = j)
            }
            
          }
          
          temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- setdiff(temp_vector_names_of_newly_discovered_nested_functions_for_inspection0, ls(all.names = TRUE))
          
          temp_current_frame <- temp_current_frame - 1
          
        }
        
        # uncollected function names are added to the global uncollected list
        temp_global_variables_uncollected <- c(temp_global_variables_uncollected, temp_vector_names_of_newly_discovered_nested_functions_for_inspection0) %>% unique
        
        # flag for further internal inspection (update to have only collected and loaded functions)
        temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- intersect(temp_vector_names_of_newly_discovered_nested_functions_for_inspection, lsf.str() %>% as.character)
        
        # collected functions added to vector_global_variables
        vector_global_variables <- c(vector_global_variables, temp_vector_names_of_newly_discovered_nested_functions_for_inspection) %>% unique
        
        temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- character()
        
      }
      
      if (length(temp_global_variables_uncollected) > 0) {
        warning(paste("Scoping has finished but there remain some uncollected variables: ", paste(temp_global_variables_uncollected, collapse = " , "), "\n"))
        
        if (.debug == TRUE) {
          message("all variables in environment right before writing to disk")
          print(ls(all.names = TRUE))
        }
        
        vector_global_variables <- setdiff(vector_global_variables, temp_global_variables_uncollected)
      }
      
      message("Writing globals to disk")
      # save(list = vector_global_variables, file = .intermediate_files_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.intermediate_files_dir, "/", i, ".qs", sep = ""))
      }
      
      # write vector_global_variables
      qs::qsave(x = vector_global_variables, file = paste(.intermediate_files_dir, "/", .job_name, "_vector_global_variables.qs", sep = ""))
      
    }
    
  }
  
  # END SCOPING ####
  
  # check if all the list elements are of equal length
  .list_length <- unique(unlist(lapply(X = .l, FUN = function(a1) {return(length(a1))} )))
  
  # print(".list_length")
  # print(.list_length)
  # print(".no_workers")
  # print(.no_workers)
  
  map_length <- .list_length
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  if (.no_workers > map_length) {
    .no_workers <- map_length
  }
  
  if (map_length == 0) {
    return(.l[[1]])
  }
  
  # no_nodes_expected_to_be_utilised <- ceiling(.no_workers/.no_workers_per_node)
  
  # preallocate worker list
  list_workers <- as.list(paste("chunk_", 1:map_length, sep = ""))
  names(list_workers) <- paste("chunk_", 1:map_length, sep = "")
  
  # CHUNKING ####
  
  if (is.null(.no_chunks) == TRUE) {
    .no_chunks <- .no_workers
  } else if (is.numeric(.no_chunks) == FALSE) {
    .no_chunks <- .no_workers
  }
  
  if (.no_chunks > map_length) {
    .no_chunks <- map_length
  }
  
  map_length <- .no_chunks
  
  list_splitting_schema <- al_splitindices(index_length = .list_length, number_of_chunks = .no_chunks)
  
  # write .l
  lapply(
    X = 1:.no_chunks, 
    FUN = function(a1) {
      
      qs::qsave(
        x = lapply(
          X = .l,
          FUN = function(b1) {
            return(
              b1[list_splitting_schema[[a1]]$start:list_splitting_schema[[a1]]$end]
            )
          } )
        , file = paste(.intermediate_files_dir, "/", .job_name, "_.l_chunk_", a1, ".qs", sep = "")
      )
      
    }
  )
  
  # write .f
  qs::qsave(x = .f, file = paste(.intermediate_files_dir, "/", .job_name, "_.f.qs", sep = ""))
  
  # END CHUNKING ####
  
  # WRITE RUNTIME SCRIPTS
  # L1: bash ... 1
  #     bash ... 2
  #     bash ... 3
  #     bash ... 4
  # L2: declare bash variables, calling Rscript <L3.R> - 1 file, L1 parameters accepted
  # L3: the actual commands in R - 1 file, parameter passthrough
  
  # function_to_run <- function(.l_current, .f, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
  
  # commandArgs(trailingOnly = TRUE)
  # [1] ..i (chunk no.)
  
  # L1 ####
  ## annoyingly, the L1 verification only works when it's called from a script, not through the nci-parallel bash loop. so fucking stupid.
  writeLines(text = paste("echo \"L1_VERIFY\" > \"", .intermediate_files_dir, "/", .job_name, "_nci_multinode_L1_verification.txt\"", sep = ""), con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L1_verification_cmd.txt", sep = ""))
  
  writeLines(text = c(paste("bash ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_L1_verification_cmd.txt", sep = ""), paste("bash ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt ", 1:.no_chunks, sep = "")), con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L1.txt", sep = ""))
  # END L1 ####
  
  # L2 ####
  write(x = paste("echo \"L2_VERIFY_INVOCATION\" > \"", .intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_invocation_chunk_\"$1\".txt\"", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""), append = FALSE)
  
  system(command = paste("printenv | grep ^PATH= | sed 's|^\\([^=]*\\)=\\(.*\\)|export \\1=\"\\2\"|g' >> ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""))
  system(command = paste("printenv | grep ^LIBRARY_PATH= | sed 's|^\\([^=]*\\)=\\(.*\\)|export \\1=\"\\2\"|g' >> ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""))
  system(command = paste("printenv | grep ^LD_LIBRARY_PATH= | sed 's|^\\([^=]*\\)=\\(.*\\)|export \\1=\"\\2\"|g' >> ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""))
  system(command = paste("printenv | grep ^R_LIBS_USER= | sed 's|^\\([^=]*\\)=\\(.*\\)|export \\1=\"\\2\"|g' >> ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""))
  system(command = paste("printenv | grep ^TMPDIR= | sed 's|^\\([^=]*\\)=\\(.*\\)|export \\1=\"\\2\"|g' >> ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""))
  
  # necessary to test if R can run before we verify
  write(x = paste("Rscript -e \"q()\"", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""), append = TRUE)
  
  write(x = paste("echo \"L2_VERIFY_RSCRIPT\" > \"", .intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_rscript_chunk_\"$1\".txt\"", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""), append = TRUE)
  
  # $1 ..i (chunk no.)
  write(x = paste("Rscript ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R $1 ", .epoch_time, " > \"", .status_messages_dir, "/", .job_name, "_chunk_\"$i\"_stdout.txt\" 2> \"", .status_messages_dir, "/", .job_name, "_chunk_\"$i\"_stderr.txt\"", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L2.txt", sep = ""), append = TRUE)
  # END L2 ####
  
  # L3 ####
  write(x = "", file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = FALSE)
  
  write(x = paste("writeLines(text = \"\", con = paste(\"", .intermediate_files_dir, "/", .job_name, "_chunk_\", commandArgs(trailingOnly = TRUE)[1], \"_exitmsg.txt\", sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste(".intermediate_files_dir <- \"", .intermediate_files_dir, "\"", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste(".vector_global_packages <- qs::qread(file = \"", .intermediate_files_dir, "/", .job_name, "_vector_global_packages.qs", "\")", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  write(x = paste(".vector_global_variables <- qs::qread(file = \"", .intermediate_files_dir, "/", .job_name, "_vector_global_variables.qs", "\")", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  write(x = paste(".l <- qs::qread(file = paste(\"", .intermediate_files_dir, "/", .job_name, "_.l_chunk_\", commandArgs(trailingOnly = TRUE)[1], \".qs\", sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  write(x = paste(".f <- qs::qread(file = \"", .intermediate_files_dir, "/", .job_name, "_.f.qs", "\")", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("lapply(.vector_global_packages, library, character.only = TRUE)", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("for (i in .vector_global_variables) {
    assign(x = i, value = qs::qread(file = paste(\"", .intermediate_files_dir, "\", i, \".qs\", sep = \"\")), envir = .GlobalEnv)
  }", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("if (", .debug, " == TRUE) {
    print(\"ls before pmap\")
    print(ls(all.names = TRUE))
  }", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("obj <- purrr::pmap(
    .l = .l,
    .f = .f,
    .progress = ", .progress, "
  )", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("qs::qsave(x = obj, file = paste(\"", .intermediate_files_dir, "/", .job_name, "_chunk_\", commandArgs(trailingOnly = TRUE)[1], \".qs\", sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("writeLines(text = \"GRACEFUL_EXIT\", con = paste(\"", .intermediate_files_dir, "/", .job_name, "_chunk_\", commandArgs(trailingOnly = TRUE)[1], \"_exitmsg.txt\", sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("system(paste(\"kill -9 \", Sys.getpid(), sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  write(x = paste("system(paste(\"kill -9 \", Sys.getpid(), \"; wait \", Sys.getpid(), \" 2>/dev/null || true; exit 0\", sep = \"\"))", sep = ""), file = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L3.R", sep = ""), append = TRUE)
  
  # END L3 ####
  
  message("Commencing computation")
  
  # INITIALISATION BEFORE LOOP
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  flag_completion <- FALSE
  
  vector_exit_codes <- NULL
  
  vector_names_of_chunks_alive <- character()
  
  vector_logical_indices_workers_completed_reported <- rep(x = FALSE, times = map_length)
  names(vector_logical_indices_workers_completed_reported) <- paste("chunk_", 1:map_length, sep = "")
  
  vector_current_chunks_spliced <- numeric()
  
  # END INITIALISATION ###
  
  system(command = paste("echo $$ > ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_mother_pid.txt; module load nci-parallel; mpirun --np ", .no_workers, " --map-by ppr:", .no_workers_per_node, ":node:PE=", .no_cores_per_worker, " --oversubscribe --bind-to none nci-parallel --dedicated --input-file ", .intermediate_files_dir, "/", .job_name, "_nci_multinode_cmdfile_L1.txt ", .epoch_time, " &", sep = ""), intern = FALSE)
  
  mpirun_mother_pid <- readLines(con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_mother_pid.txt", sep = ""))
  
  # formulate our own exit codes
  origin_hostname_raw <- system(command = "hostname", intern = TRUE)
  
  ## retrieve system process status
  ## KNOWN BUG: PBS_JOBID is old if you don't quit session before running
  pbs_jobname <- system(command = "echo $PBS_JOBID", intern = TRUE)
  
  vector_node_hostnames <- system(command = paste("qstat -an1 ", pbs_jobname, " | awk 'NR == 6{print $12}' | awk '{gsub(/\\/.\\*[0-9]+\\+{0,1}/, \" \"); print $0}'", sep = ""), intern = TRUE)
  vector_node_hostnames <- unlist(strsplit(trimws(vector_node_hostnames), "\\s"), recursive =  TRUE)
  message(paste("Running job on hosts: ", paste(vector_node_hostnames, collapse = ","), sep = ""))
  
  origin_hostname_short <- vector_node_hostnames[unlist(lapply(X = vector_node_hostnames, FUN = function(a1) {return(grep(x = origin_hostname_raw, pattern = a1))}))]
  global_origin_hostname_short <<- origin_hostname_short
  
  # this will hang if L1 cant verify (if nci-parallel is not launching the L1 script)
  message("Awaiting L1 verification")
  # we can wait for a bit for the verification to fully write if disk writes are very slow
  # after a while, we lose patience
  temp_patience_counter <- 0
  while (!file.exists(paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L1_verification.txt", sep = ""))) {
  }
  while (readLines(con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L1_verification.txt", sep = "")) != "L1_VERIFY" & temp_patience_counter < .patience) {
    Sys.sleep(1)
    temp_patience_counter <- temp_patience_counter + 1
  }
  message("L1 verification complete")
  
  # construct the master table that records process status
  df_master_process_status <- data.frame(
    "chunkname" = names(list_workers),
    "pid" = NA,
    "hostname" = NA,
    "exitmsg" = NA,
    "verify_invocation" = NA,
    "verify_rscript" = NA,
    "spawned" = FALSE,
    "current_process_status_code" = NA,
    "patience" = 0,
    "exit_code" = -2
  )
  
  logical_r_processes_have_spawned <- FALSE
  
  while(all(vector_exit_codes <= 0)) {
    
    # we must read exitmsgs BEFORE polling processes because exitmsgs come after process launch
    vector_exitmsg <- unlist(lapply(
      X = names(list_workers),
      FUN = function(a1) {
        
        if (file.exists(paste(.intermediate_files_dir, "/", .job_name, "_", a1, "_exitmsg.txt", sep = ""))) {
          return(readLines(con = paste(.intermediate_files_dir, "/", .job_name, "_", a1, "_exitmsg.txt", sep = "")))
        } else {
          return("missing")
        }
        
      } ))
    
    # add exitmsg to master process table
    df_master_process_status$exitmsg <- vector_exitmsg
    
    vector_multinode_L2_verify_invocation <- unlist(lapply(
      X = names(list_workers),
      FUN = function(a1) {
        
        if (file.exists(paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_invocation_", a1, ".txt", sep = ""))) {
          return(readLines(con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_invocation_", a1, ".txt", sep = "")))
        } else {
          return("missing")
        }
        
      } ))
    
    vector_multinode_L2_verify_rscript <- unlist(lapply(
      X = names(list_workers),
      FUN = function(a1) {
        
        if (file.exists(paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_rscript_", a1, ".txt", sep = ""))) {
          return(readLines(con = paste(.intermediate_files_dir, "/", .job_name, "_nci_multinode_L2_verification_rscript_", a1, ".txt", sep = "")))
        } else {
          return("missing")
        }
        
      } ))
    
    # WRANGLE LIVE PROCESS DATA
    ## this a snapshot of all processes running on the nodes
    list_df_ps <- lapply(
      X = vector_node_hostnames, 
      FUN = function(a1) {
        
        # DEBUG ###
        # a1 <- vector_node_hostnames[1]
        ###########
        
        # we cannot use ssh when checking on the same node
        if (a1 == origin_hostname_short) {
          vector_ps_status <- trimws(system(command = "ps -eo pid,ppid,pgid,s,cmd", intern = TRUE))[-1]
        } else {
          vector_ps_status <- trimws(system(command = paste("ssh -o StrictHostKeyChecking=no ", a1, " 'ps -eo pid,ppid,pgid,s,cmd'", sep = ""), intern = TRUE))[-1]
        }
        
        # we are doing this because there's no way to split by spaces without also separating the CMD column. no, tab separation attempts don't work.
        L1_df_ps <- data.frame(
          "pid" = gsub(x = vector_ps_status, pattern = "^([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+(.+)$", replacement = "\\1"),
          "ppid" = gsub(x = vector_ps_status, pattern = "^([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+(.+)$", replacement = "\\2"),
          "pgid" = gsub(x = vector_ps_status, pattern = "^([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+(.+)$", replacement = "\\3"),
          "s" = gsub(x = vector_ps_status, pattern = "^([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+(.+)$", replacement = "\\4"),
          "cmd" = gsub(x = vector_ps_status, pattern = "^([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+([^ ]+)\\s+(.+)$", replacement = "\\5"),
          "hostname" = a1
        )
        
        return(L1_df_ps)
        
      } )
    
    # now, our pidtable has an additional hostname column
    # df_ps_child_processes is a table of all child R/rsession processes with pid, status and hostname
    ## make a base version of dplyr::bind_rows
    df_ps <- list_df_ps[[1]][0, ]
    for (i in (1:length(list_df_ps))) {
      df_ps <- rbind(df_ps, list_df_ps[[i]])
    }
    
    # simpler way to identify child processes
    df_ps_child_processes <- df_ps[grep(x = df_ps$cmd, pattern = paste(".*_nci_multinode_cmdfile_L3.R --args ([^ ]+) ", .epoch_time, sep = "")), ]
    
    df_master_process_status$current_process_status_code <- NA
    
    # infer chunk number from the chunking command
    if (nrow(df_ps_child_processes) > 0) {
      df_ps_child_processes$chunkname <- gsub(x = df_ps_child_processes$cmd, pattern = paste(".*_nci_multinode_cmdfile_L3.R --args ([^ ]+) ", .epoch_time, sep = ""), replacement = "chunk_\\1")
      
      for (temp_chunkname in df_master_process_status$chunkname[df_master_process_status$chunkname %in% df_ps_child_processes$chunkname]) {
        df_master_process_status[df_master_process_status$chunkname == temp_chunkname, "pid"] <- df_ps_child_processes[df_ps_child_processes$chunkname == temp_chunkname, ]$pid
        df_master_process_status[df_master_process_status$chunkname == temp_chunkname, "hostname"] <- df_ps_child_processes[df_ps_child_processes$chunkname == temp_chunkname, ]$hostname
        df_master_process_status[df_master_process_status$chunkname == temp_chunkname, "current_process_status_code"] <- df_ps_child_processes[df_ps_child_processes$chunkname == temp_chunkname, ]$s
        df_master_process_status[df_master_process_status$chunkname == temp_chunkname, "spawned"] <- TRUE
      }
      
    }
    
    # make df_process_status chunknames consistent with names(list_workers)
    # df_process_status <- merge(x = data.frame("chunkname" = names(list_workers)), y = df_ps_child_processes, by = "chunkname", all.x = TRUE, sort = FALSE)
    # df_process_status <- df_process_status[match(intersect(names(list_workers), df_process_status$chunkname), df_process_status$chunkname), ]
    # df_process_status[is.na(df_process_status$s), "s"] <- "missing"
    # df_process_status$isalive <- !grepl(x = df_process_status$s, pattern = "Z|X|T|t|missing", ignore.case = FALSE)
    
    # END WRANGLE LIVE PROCESS DATA ###
    
    ## add in process information as it becomes available
    ### current_process_status_code is a dynamic column
    
    # update exit codes
    # patience deals with hangs before the rscript call
    list_updated_exit_codes_and_patience <- mapply(
      "L1_chunkname" = df_master_process_status$chunkname,
      "L1_pid" = df_master_process_status$pid,
      "L1_hostname" = df_master_process_status$hostname,
      "L1_verify_invocation" = df_master_process_status$verify_invocation,
      "L1_verify_rscript" = df_master_process_status$verify_rscript,
      "L1_exitmsg" = df_master_process_status$exitmsg,
      "L1_spawned" = df_master_process_status$spawned,
      "L1_current_process_status_code" = df_master_process_status$current_process_status_code,
      "L1_patience" = df_master_process_status$patience,
      SIMPLIFY = FALSE,
      FUN = function(L1_chunkname, L1_pid, L1_hostname, L1_verify_invocation, L1_verify_rscript, L1_exitmsg, L1_spawned, L1_current_process_status_code, L1_patience) {
        
        if (L1_exitmsg == "missing" & L1_spawned == FALSE) {
          return(
            list(
              "L1_exit_code" = -2,
              "patience" = 0
            )
          )
        } else if (L1_exitmsg == "" & L1_spawned == TRUE & !grepl(x = L1_current_process_status_code, pattern = "Z|X|T|t", ignore.case = FALSE)) {
          return(
            list(
              "L1_exit_code" = -1,
              "patience" = 0
            )
          )
        } else if (L1_exitmsg == "GRACEFUL_EXIT") {
          return(
            list(
              "L1_exit_code" = 0,
              "patience" = 0
            )
          )
        } else if (L1_exitmsg == "" & is.na(L1_current_process_status_code) == TRUE) {
          return(
            list(
              "L1_exit_code" = 1,
              "patience" = 0
            )
          )
        } else {
          # the remaining cases are vast majority due to L1_exitmsg == "missing" & L1_spawned == TRUE. that is a patience issue
          return(
            list(
              "L1_exit_code" = -2,
              "patience" = L1_patience + 1
            )
          )
        }
        
      } )
    
    global_list_updated_exit_codes_and_patience <<- list_updated_exit_codes_and_patience
    
    vector_updated_exit_codes <- unlist(lapply(X = list_updated_exit_codes_and_patience, FUN = function(a1) {return(a1$L1_exit_code)}))
    vector_updated_patience <- unlist(lapply(X = list_updated_exit_codes_and_patience, FUN = function(a1) {return(a1$patience)}))
    
    # write new params
    global_df_master_process_status <<- df_master_process_status
    
    df_master_process_status$exit_code <- vector_updated_exit_codes
    df_master_process_status$patience <- vector_updated_patience
    
    ## deal with the chunks we expect to be still running or not
    if (length(vector_current_chunks_spliced) > 0) {
      vector_names_of_current_chunks_spliced <- paste("chunk_", vector_current_chunks_spliced, sep = "")
      
      if (.debug == TRUE) {
        global_vector_current_chunks_spliced <<- vector_current_chunks_spliced
      }
      
    } else {
      vector_names_of_current_chunks_spliced <- character()
    }
    
    vector_names_of_chunks_alive <- df_master_process_status[df_master_process_status$exit_code == -1,]$chunkname
    
    vector_exit_codes <- df_master_process_status$exit_code
    names(vector_exit_codes) <- names(list_workers)
    
    if (.debug == TRUE) {
      print("df_process_status")
      print(df_process_status)
      
      global_df_process_status <<- df_process_status
      global_vector_names_of_chunks_alive <<- vector_names_of_chunks_alive
      global_vector_names_of_current_chunks_spliced <<- vector_names_of_current_chunks_spliced
      global_vector_exit_codes <<- vector_exit_codes
      global_list_workers <<- list_workers
    }
    
    # detect error
    ## if any error, stop all
    if (any(vector_exit_codes > 0) | any(df_master_process_status$patience == .patience)) {
      print("Process status per chunk")
      print(vector_exit_codes)
      # for mpirun, killing the mother pid will kill all the spawns in OTHER nodes. this saves us having to ssh into each node and kill.
      system(command = paste("kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste("mpirun.+nci_multinode_cmdfile_L1.txt ", .epoch_time, sep = "")),]$pid, collapse = " "), "; kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste(".*_nci_multinode_cmdfile_L3.R --args ([^ ]+) ", .epoch_time, sep = "")),]$pid, collapse = " "), sep = ""), ignore.stdout = TRUE, ignore.stderr = TRUE)
      options(warning.length = 8170)
      if (any(vector_exit_codes > 0)) {
        stop(paste("Exit status failure received on chunks: ", paste(gsub(x = names(vector_exit_codes[vector_exit_codes > 0]), pattern = "chunk_", replacement = ""), collapse = ","), ". \n\nPlease run the following command in order to view the outputs of each worker: \n", "for i in ", paste(names(vector_exit_codes[vector_exit_codes > 0]), collapse = " "), "; do tail -n 20 ", "\"", .status_messages_dir, "\"", "$i\"_stdout.txt\" ", "\"", .status_messages_dir, "\"", "$i\"_stderr.txt\"; done", sep = ""))
      } else if (any(df_master_process_status$patience == .patience)) {
        stop(paste("Patience exhausted on chunks: ", paste(gsub(x = df_master_process_status[df_master_process_status$patience == .patience, ]$chunkname, pattern = "chunk_", replacement = ""), collapse = ","), ". \n\nPlease run the following command in order to view the outputs of each worker: \n", "for i in ", paste(df_master_process_status[df_master_process_status$patience == .patience, ]$chunkname, collapse = " "), "; do tail -n 20 ", "\"", .status_messages_dir, "\"", "$i\"_stdout.txt\" ", "\"", .status_messages_dir, "\"", "$i\"_stderr.txt\"; done", sep = ""))
      }
    }
    
    vector_logical_indices_workers_completed_reported <- vector_exit_codes == 0
    
    if (.debug == TRUE) {
      print("vector_logical_indices_workers_completed_reported")
      print(vector_logical_indices_workers_completed_reported)
      
      global_vector_logical_indices_workers_completed_reported <<- vector_logical_indices_workers_completed_reported
    }
    
    if (flag_completion == TRUE) {
      break()
    } else if (flag_completion == FALSE) {
      
      # SPLICER
      
      # initialise list_result if it hasnt already been created
      if (length(ls(pattern = "^list_result$")) == 0) {
        list_result <- list()
      }
      
      if (length(vector_current_chunks_spliced) < map_length) {
        
        vector_completed_chunks <- which(vector_logical_indices_workers_completed_reported == TRUE)
        
        if (.debug == TRUE) {
          print("vector_completed_chunks")
          print(vector_completed_chunks)
          global_vector_completed_chunks <<- vector_completed_chunks
        }
        
        vector_chunks_to_be_spliced <- sort(unique(setdiff(vector_completed_chunks, vector_current_chunks_spliced)))
        
        if (.debug == TRUE) {
          print("vector_chunks_to_be_spliced")
          print(vector_chunks_to_be_spliced)
          global_vector_chunks_to_be_spliced <<- vector_chunks_to_be_spliced
        }
        
        # splice terms in order
        if (length(vector_chunks_to_be_spliced) > 0) {
          
          for (.j in vector_chunks_to_be_spliced) {
            
            list_result[[.j]] <- qs::qread(file = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", .j, ".qs", sep = ""))
            
          }
          
          vector_current_chunks_spliced <- sort(unique(c(vector_current_chunks_spliced, vector_chunks_to_be_spliced)))
          
        }
        
        if (.debug == TRUE) {
          Sys.sleep(1)
        }
        
      }
      
      # END SPLICER ###
      
      # deal with completion flag
      flag_completion <- (length(vector_current_chunks_spliced) == map_length) & (length(which(vector_logical_indices_workers_completed_reported)) == map_length)
      
    }
    
    # an actually useful progress bar although rudimentary
    if (length(ls(pattern = "^list_result$")) == 0) {
      cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%)", sep = ""))
    } else {
      cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%); Splice progress: ", length(vector_current_chunks_spliced), "/", map_length, sep = ""))
    }
    
    if (.debug == TRUE) {
      Sys.sleep(1)
    } else {
      Sys.sleep(5)
    }
    
  }
  
  # detect error
  ## if any error, stop all
  if (any(vector_exit_codes > 0)) {
    print("Process status per chunk")
    print(vector_exit_codes)
    system(command = paste("kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste("mpirun.+nci_multinode_cmdfile_L1.txt ", .epoch_time, sep = "")),]$pid, collapse = " "), "; kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste(".*_nci_multinode_cmdfile_L3.R --args ([^ ]+) ", .epoch_time, sep = "")),]$pid, collapse = " "), sep = ""), ignore.stdout = TRUE, ignore.stderr = TRUE)
    # spew out the last 20 status lines for easy debugging
    # lapply(X = names(vector_exit_codes)[vector_exit_codes > 0], FUN = function(a1) {warning(paste(a1, " stdout (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stdout.txt", sep = ""), sep = ""))); warning(paste(a1, " stderr (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stderr.txt", sep = ""), sep = "")))} )
    options(warning.length = 8170)
    stop(paste("Exit status failure received on chunks: ", paste(gsub(x = names(vector_exit_codes[vector_exit_codes > 0]), pattern = "chunk_", replacement = ""), collapse = ","), ". \n\nPlease run the following command in order to view the outputs of each worker: \n", "for i in ", paste(names(vector_exit_codes[vector_exit_codes > 0]), collapse = " "), "; do tail -n 20 ", "\"", .status_messages_dir, "\"", "$i\"_stdout.txt\" ", "\"", .status_messages_dir, "\"", "$i\"_stderr.txt\"; done", sep = ""))
    # } else {
    # gotta leave the process hanging for a bit so we can verify it finished correctly. THEN we manually kill the child.
    # lapply(X = list_workers[names(vector_exit_codes)[vector_exit_codes == 0]], FUN = function(a1) {system(command = paste("kill -9 ", mpirun_mother_pid, sep = ""))} )
  }
  
  if (.keep_intermediate_files == FALSE) {
    # unlink(list.files(path = .intermediate_files_dir, pattern = paste(.job_name, "_chunk_.*.rds", sep = ""), full.names = TRUE ), recursive = TRUE)
    # unlink(list.files(path = .temp_dir, pattern = paste(.job_name, ".rdata", sep = ""), full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .intermediate_files_dir, pattern = ".*", full.names = TRUE ), recursive = TRUE)
  }
  
  system(command = paste("kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste("mpirun.+nci_multinode_cmdfile_L1.txt ", .epoch_time, sep = "")),]$pid, collapse = " "), "; kill -9 ", paste(df_ps[(df_ps$hostname == origin_hostname_short) & grepl(x = df_ps$cmd, pattern = paste(".*_nci_multinode_cmdfile_L3.R --args ([^ ]+) ", .epoch_time, sep = "")),]$pid, collapse = " "), sep = ""), ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  cat("\n")
  
  # unchunkify
  list_result <- unlist(list_result, recursive = FALSE, use.names = TRUE)
  
  return(list_result)
  
}
