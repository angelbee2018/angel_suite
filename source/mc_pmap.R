mc_pmap <- function(
    .l, .f, 
    .no_workers = 1, .no_chunks = 1, .splicing_order = "ordered", 
    .job_name = NULL, 
    .globals_save_dir = NULL, .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL, 
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE, 
    .status_messages_dir = NULL, .progress = TRUE,
    .debug = FALSE) {
  
  # DEBUG ###
  # .l = list(
  #   "b1" = 1:20
  # )
  # .no_workers = 20
  # .no_chunks <- 20
  # .globals_mode = "user"
  # .re_export = TRUE
  # .globals_save_dir = paste(tempdir(), "/tempdata.rdata", sep = "")
  # .intermediate_files_dir = tempdir()
  # .user_global_objects = c()
  # .status_messages_dir = paste(tempdir(), sep = "")
  # .progress <- TRUE
  # .job_name = "test"
  # .splicing_order = "ordered"
  # .f = function(b1) {set.seed(b1);Sys.sleep(runif(n = 1, min = 3, max = 5)); return(list(LETTERS[b1]))}
  # .debug <- FALSE
  ###########

  for (i in c("parallel", "utils", "lubridate", "qs")) { 
    
    if (require(i, character.only = TRUE) == FALSE) {
      stop(paste("Package \"", i, "\" not found. Please install using `install.packages` or `BiocManager::install`", sep = ""))
    }
    
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    if (type.convert(system("ulimit -n", intern = TRUE), as.is = TRUE) < 65536) {
      warning(paste("System max. open file limit is less than the recommended 65536. Current limit is set to: ", system("ulimit -n", intern = TRUE), ". To fix this, please re-run R from bash terminal after having set `ulimit -n 65536` to avoid possible errors with large jobs and/or large number of workers/chunks.", sep = ""))
    }
  }
  
  .epoch_time <- as.numeric(Sys.time())*1E5
  .temp_dir <- paste(tempdir(), "_", .epoch_time, "/", sep = "")
  
  if (is.null(.job_name)) {
    .job_name <- paste("round_robin_pmap_callr_", .epoch_time, sep = "")
  }
  
  message(paste("Job name: ", .job_name, sep = ""))
  
  if (is.null(.globals_save_dir)) {
    .globals_save_dir <- .temp_dir
  }
  if (!dir.exists(.globals_save_dir)) {
    dir.create(.globals_save_dir, recursive = TRUE)
  }
  
  message(paste("Globals will be saved to disk at: ", .globals_save_dir, sep = ""))
  
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
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(lapply(X = .l, FUN = function(a1) {return(length(a1))} )))
  
  if (.no_workers > map_length) {
    .no_workers <- map_length
  }
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  if (map_length == 0) {
    return(.l[[1]])
  }
  
  # preallocate worker list
  list_workers <- list()
  
  # CHUNKING ###
  
  if (is.null(.no_chunks) == TRUE) {
    .no_chunks <- .no_workers
  } else if (is.numeric(.no_chunks) == FALSE) {
    .no_chunks <- .no_workers
  }
  
  if (.no_chunks > map_length) {
    .no_chunks <- map_length
  }
  
  map_length <- .no_chunks
  
  # END CHUNKING ###

# define function to be fun
          function_to_run <- function(.l_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, ..i, .status_messages_dir_stdout, .status_messages_dir_stderr) {
            
            writeLines(text = "", con = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", ..i, "_exitmsg.txt", sep = ""))
            
            .status_messages_dir_stdout_con <- file(.status_messages_dir_stdout, open = "a")
            sink(file = .status_messages_dir_stdout_con, split = TRUE, append = FALSE, type = "output")
            
            .status_messages_dir_stderr_con <- file(.status_messages_dir_stderr, open = "a")
            sink(file = .status_messages_dir_stderr_con, split = FALSE, append = FALSE, type = "message")
            
            # lapply(vector_global_packages, library, character.only = TRUE)
            # for (i in vector_global_variables) {
            #   assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
            # }
            
            # load(file = .globals_save_dir, envir = .GlobalEnv)
            
            if (.debug == TRUE) {
              print("ls before pmap")
              print(ls())
            }
            
            obj <- purrr::pmap(
              .l = .l_current,
              .f = .f,
              .progress = .progress
            )
            
            # saveRDS(object = obj, file = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", ..i, ".rds", sep = ""))
            qs::qsave(x = obj, file = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", ..i, ".qs", sep = ""))
            
            writeLines(text = "GRACEFUL EXIT", con = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", ..i, "_exitmsg.txt", sep = ""))
            
            system(paste("kill -9 ", Sys.getpid(), sep = ""))
            
          }
    
  message("Commencing computation")
      
  # keep track of iteration progress
  ## initialise
  current_map_index <- 0
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  flag_completion <- FALSE
  
  vector_exit_codes <- NULL
  
  vector_logical_indices_workers_completed_reported <- rep(x = FALSE, times = map_length)
  names(vector_logical_indices_workers_completed_reported) <- paste("chunk_", 1:map_length, sep = "")
  
  vector_current_chunks_spliced <- numeric()
  
  while(all(vector_exit_codes <= 0)) {
    
    # formulate our own exit codes

    if (current_map_index > 0) {

        ## retrieve system process status
    vec_topresult <- trimws(system(command = "ps -eo pid,s", intern =  TRUE))
    df_topresult <- as.data.frame(t(as.data.frame(strsplit(vec_topresult[2:length(vec_topresult)], split = "\\s+"))))
    colnames(df_topresult) <- unlist(strsplit(vec_topresult[1], split = "\\s+"))
    
    ## retrieve our list of assigned pids per chunk
    df_pidtable <- data.frame("chunkname" = names(list_workers), "PID" = unlist(lapply(X = list_workers, FUN = function(a1) {a1$pid})))

    ## deal with the chunks we expect to be still running or not
    if (length(vector_current_chunks_spliced) > 0) {
      vector_names_of_current_chunks_spliced <- paste("chunk_", vector_current_chunks_spliced, sep = "")
      
      if (.debug == TRUE) {
        global_vector_current_chunks_spliced <<- vector_current_chunks_spliced
      }
      
    } else {
      vector_names_of_current_chunks_spliced <- character()
    }
    
    # tibble_process_status <- dplyr::left_join(df_pidtable, df_topresult, by = "PID")
    df_process_status <- merge(x = df_pidtable, y = df_topresult, by = "PID", all.x = TRUE)
    
    df_process_status[is.na(df_process_status$S), "S"] <- "missing"
    df_process_status$isalive <- !grepl(x = df_process_status$S, pattern = "Z|X|T|t|missing", ignore.case = FALSE)
    
    vector_chunks_alive <- df_process_status[df_process_status$isalive == TRUE, ]$chunkname
    
    vector_exit_codes <- mapply(
      "a1" = list_workers,
      "a2" = names(list_workers),
      "a3" = df_process_status$isalive,
      FUN = function(a1, a2, a3) {
        
        # logical_process_is_alive <- type.convert(system(command = paste("ps -r ", a1$pid, " | wc -l", sep = ""), intern = TRUE), as.is = TRUE) > 1
        
        logical_process_is_alive <- a3
        
        exitmsg_path <- paste(.intermediate_files_dir, "/", .job_name, "_", a2, "_exitmsg.txt", sep = "")
        
        if (file.exists(exitmsg_path)) {
          exitmsg <- readLines(con = exitmsg_path)
          
          if (exitmsg == "GRACEFUL EXIT" | a2 %in% vector_names_of_current_chunks_spliced) {
            return(0)
          } else if (logical_process_is_alive == TRUE) {
            return(-1)
          } else if (logical_process_is_alive == FALSE) {
            return(1)
          }
          
        } else {
          return("2")
        }
        
      } )
    
    names(vector_exit_codes) <- names(list_workers)
        
    } else {
      vector_exit_codes <- numeric()
    }   
    
    if (.debug == TRUE) {
      print("df_process_status")
      print(df_process_status)
      
      global_df_process_status <<- df_process_status
      global_vector_chunks_alive <<- vector_chunks_alive
      global_vector_names_of_current_chunks_spliced <<- vector_names_of_current_chunks_spliced
      global_vector_exit_codes <<- vector_exit_codes
      global_list_workers <<- list_workers
    }
    
    # detect error
    ## if any error, stop all
    if (any(vector_exit_codes > 0)) {
      print("Process status per chunk")
      print(vector_exit_codes)
      lapply(X = list_workers, FUN = function(a1) {system(paste("kill -9 ", a1$pid, sep = ""))} )
      # spew out the last 20 status lines for easy debugging
      # lapply(X = names(vector_exit_codes)[vector_exit_codes > 0], FUN = function(a1) {warning(paste(a1, " stdout (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stdout.txt", sep = ""), sep = ""))); warning(paste(a1, " stderr (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stderr.txt", sep = ""), sep = "")))} )
      options(warning.length = 8170)
      stop(paste("Exit status failure received on chunks: ", paste(gsub(x = names(vector_exit_codes[vector_exit_codes > 0]), pattern = "chunk_", replacement = ""), collapse = ","), ". \n\nPlease run the following command in order to view the outputs of each worker: \n", "for i in ", paste(names(vector_exit_codes[vector_exit_codes > 0]), collapse = " "), " ; do tail -n 20 ", "\"", .status_messages_dir, "\"", "/$i\"_stdout.txt\" ", "\"", .status_messages_dir, "\"", "/$i\"_stderr.txt\"; done", sep = ""))
      # } else {
      # gotta leave the process hanging for a bit so we can verify it finished correctly. THEN we manually kill the child.
      # lapply(X = list_workers[names(vector_exit_codes)[vector_exit_codes == 0]], FUN = function(a1) {system(paste("kill -9 ", a1$pid, sep = ""))} )
    }
    
    vector_logical_indices_workers_completed_reported <- vector_exit_codes == 0
    
    if (.debug == TRUE) {
      print("vector_logical_indices_workers_completed_reported")
      print(vector_logical_indices_workers_completed_reported)
      
      global_vector_logical_indices_workers_completed_reported <<- vector_logical_indices_workers_completed_reported
    }
    
    # round robin action until the list is fully mapped
    if (flag_completion == TRUE) {
      break()
    } else if (flag_completion == FALSE) {
      
      # add more processes as long as we need to fill up worker slots and we're not at the end of the list yet
      number_of_processes_alive <- length(vector_chunks_alive)
      if (number_of_processes_alive < .no_workers & current_map_index < map_length) {
        
        new_map_start <- current_map_index + 1
        new_map_end <- min(c(current_map_index + .no_workers - number_of_processes_alive, map_length))
        
        for (..i in (new_map_start):min(c(new_map_end, map_length))) {
                    
          # define chunk for the target of mapping operation
          .l_current <- lapply(X = .l, FUN = function(a1) {return(a1[parallel::splitIndices(nx = length(.l[[1]]), ncl = .no_chunks)[[..i]]])} )
          
          .status_messages_dir_stdout <- paste(.status_messages_dir, "/chunk_", ..i, "_stdout.txt", sep = "")
          .status_messages_dir_stderr <- paste(.status_messages_dir, "/chunk_", ..i, "_stderr.txt", sep = "")
          
          list_workers[[..i]] <- parallel:::mcparallel(
              expr = function_to_run(".l_current" = .l_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "..i" = ..i, ".status_messages_dir_stdout" = .status_messages_dir_stdout, ".status_messages_dir_stderr" = .status_messages_dir_stderr), 
              detached = TRUE
            )
                    
          names(list_workers)[..i] <- paste("chunk_", ..i, sep = "")
          
        }

          current_map_index <- length(list_workers)
        
      }
      
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
        
        # consider order
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
      cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%)", sep = ""))
    } else {
      cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%); Splice progress: ", length(vector_current_chunks_spliced), "/", map_length, sep = ""))
    }
    
    if (.debug == TRUE) {
      Sys.sleep(1)
    } else {
      Sys.sleep(0.1)
    }
    
  }
  
  # detect error
  ## if any error, stop all
  if (any(vector_exit_codes > 0)) {
    print("Process status per chunk")
    print(vector_exit_codes)
    lapply(X = list_workers, FUN = function(a1) {system(paste("kill -9 ", a1$pid, sep = ""))} )
    # spew out the last 20 status lines for easy debugging
    # lapply(X = names(vector_exit_codes)[vector_exit_codes > 0], FUN = function(a1) {warning(paste(a1, " stdout (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stdout.txt", sep = ""), sep = ""))); warning(paste(a1, " stderr (last 20 lines)", sep = "")); warning(system(command = paste("tail -n 20 ", paste(.status_messages_dir, "/", a1, "_stderr.txt", sep = ""), sep = "")))} )
    options(warning.length = 8170)
    stop(paste("Exit status failure received on chunks: ", paste(gsub(x = names(vector_exit_codes[vector_exit_codes > 0]), pattern = "chunk_", replacement = ""), collapse = ","), ". \n\nPlease run the following command in order to view the outputs of each worker: \n", "for i in ", paste(names(vector_exit_codes[vector_exit_codes > 0]), collapse = " "), " ; do tail -n 20 ", "\"", .status_messages_dir, "\"", "/$i\"_stdout.txt\" ", "\"", .status_messages_dir, "\"", "/$i\"_stderr.txt\"; done", sep = ""))
    # } else {
    # gotta leave the process hanging for a bit so we can verify it finished correctly. THEN we manually kill the child.
    # lapply(X = list_workers[names(vector_exit_codes)[vector_exit_codes == 0]], FUN = function(a1) {system(paste("kill -9 ", a1$pid, sep = ""))} )
  }
  
  if (.keep_intermediate_files == FALSE) {
    # unlink(list.files(path = .intermediate_files_dir, pattern = paste(.job_name, "_chunk_.*.rds", sep = ""), full.names = TRUE ), recursive = TRUE)
    # unlink(list.files(path = .temp_dir, pattern = paste(.job_name, ".rdata", sep = ""), full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .globals_save_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .intermediate_files_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
  }
  
  lapply(X = list_workers, FUN = function(a1) {system(paste("kill -9 ", a1$pid, sep = ""))} )
  rm(list = ls(pattern = "chunk_"))
  rm(list_workers)
  
  cat("\n")
  
  # unchunkify
  list_result <- unlist(list_result, recursive = FALSE, use.names = TRUE)
  
  return(list_result)
  
}
