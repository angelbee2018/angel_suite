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
  # .f = function(b1) {set.seed(b1);Sys.sleep(runif(n = 1, min = 10, max = 15)); return(list(LETTERS[b1]))}
  ###########
  
  for (i in c("globals", "callr", "purrr", "parallel", "magrittr", "utils", "tibble", "dplyr", "lubridate", "qs")) { 
    
    if (require(i, character.only = TRUE) == FALSE) {
      stop(paste("Package \"", i, "\" not found. Please install using `install.packages` or `BiocManager::install`", sep = ""))
    }
    
  }
  
  if (Sys.info()["sysname"] == "Linux") {
    if (system("ulimit -n", intern = TRUE) %>% type.convert(as.is = TRUE) < 65536) {
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
  
  # SCOPING ###
  
  vector_global_packages <- globals::packagesOf(globals::globalsOf(.f, mustExist = FALSE)) %>% setdiff(., c("base", "rlang"))
  
  if (.debug == TRUE) {
    print("vector_global_packages at the surface level of the function")
    print(vector_global_packages)
  }
  
  # save the current env into temp path
  if (.re_export == FALSE & dir.exists(.globals_save_dir)) {
    
    message("Temp object data file found. Will load that instead of re-exporting")
    
  } else {
    
    if (.globals_mode == "global") {
      message("Writing globals to disk")
      save.image(file = .globals_save_dir, compress = .globals_save_compress)
    } else if (.globals_mode == "auto") {
      
      vector_global_variables <- c(globals::findGlobals(.f), .user_global_objects)
      
      if (.debug == TRUE) {
        print(vector_global_variables)
      }
      
      temp_global_variables_uncollected <- setdiff(vector_global_variables, ls())
      
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
        
        temp_global_variables_uncollected <- setdiff(temp_global_variables_uncollected, ls())
        
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
          
          temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- setdiff(temp_vector_names_of_newly_discovered_nested_functions_for_inspection0, ls())
          
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
          print(ls())
        }
        
        vector_global_variables <- setdiff(vector_global_variables, temp_global_variables_uncollected)
      }
      
      message("Writing globals to disk")
      # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.globals_save_dir, i, ".qs", sep = ""))
      }
      
    }
    
  }
  
  # END SCOPING ###
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(purrr::map(.x = .l, .f = ~length(.x))))
  
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
  
  message("Commencing computation")
  
  # initial allocation of tasks to maximum number of workers
  for (..i in 1:.no_workers) {
    
    # define chunk for the target of mapping operation
    .l_current <- purrr::map(.x = .l, .f = ~.x[parallel::splitIndices(nx = length(.l[[1]]), ncl = .no_chunks)[[..i]]])
    
    # define function to be fun
    function_to_run <- function(.l_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
      
      lapply(vector_global_packages, library, character.only = TRUE)
      
      for (i in vector_global_variables) {
        assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
      }
      
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
      
      return("GRACEFUL EXIT")
      
    }
    
    # assign background worker
    assign(
      x = paste("chunk_", ..i, sep = ""), 
      value = parallel:::mcparallel(
        expr = function_to_run(".l_current" = .l_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i)
      )
    )
    
    parallel:::mccollect()
    
    list_workers <- purrr::splice(list_workers, get(paste("chunk_", ..i, sep = "")))
    names(list_workers)[length(list_workers)] <- paste("chunk_", ..i, sep = "")
    
  }
  
  global_list_workers_initial <<- list_workers
  
  # keep track of iteration progress
  current_map_index <- .no_workers
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  flag_completion <- FALSE
  
  vector_exit_statuses <- 0
  
  while(all(vector_exit_statuses == 0)) {
    
    vector_logical_indices_workers_completed_reported <- unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0
    
    if (.debug == TRUE) {
      print("vector_logical_indices_workers_completed_reported")
      print(vector_logical_indices_workers_completed_reported)
      
      global_vector_logical_indices_workers_completed_reported <<- vector_logical_indices_workers_completed_reported
    }
    
    # check on process status and write progress file
    list_logical_any_error <- purrr::map2(
      .x = list_workers, 
      .y = names(list_workers), 
      .f = function(a1, a2) {
        
        # DEBUG ###
        # a1 <- list_workers[[8]]
        # a2 <- names(list_workers) %>% .[[8]]
        ###########
        
        stdout_lines <- a1$read_output_lines()
        stderr_lines <- a1$read_error_lines()
        
        if (length(stdout_lines) != 0) {
          write(x = paste(Sys.time()), file = paste(.status_messages_dir, .job_name, "_", a2, "_stdout.txt", sep = ""), append = TRUE)
          write.table(x = stdout_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stdout.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE) %>% suppressWarnings()
        }
        
        if (length(stderr_lines) != 0) {
          write(x = paste(Sys.time()), file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE)
          write.table(x = stderr_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE) %>% suppressWarnings()
        }
        
        return(
          grepl(x = stderr_lines, pattern = "^Error in", ignore.case = FALSE)
        )
        
      } )
    
    # check for errors that don't show up in process exit status
    if (any(unlist(list_logical_any_error)) == TRUE) {
      purrr::map(.x = list_workers, .f = ~.x$signal(9))
      stop(print(paste("Stop error received on chunks:", paste(names(list_workers)[which(unlist(list_logical_any_error))], collapse = ", "))))
    }
    
    # round robin action until the list is fully mapped
    if (flag_completion == TRUE) {
      break()
    } else if (flag_completion == FALSE) {
      
      number_of_processes_alive <- length(which(unlist(purrr::map(.x = list_workers, .f = ~.x$is_alive()))))
      
      # add more processes as long as we need to fill up worker slots and we're not at the end of the list yet
      if (number_of_processes_alive < .no_workers & current_map_index < map_length) {
        
        new_map_start <- current_map_index + 1
        new_map_end <- min(c(current_map_index + .no_workers - number_of_processes_alive, map_length))
        
        current_map_index <- new_map_end
        
        for (..i in (new_map_start):min(c(new_map_end, map_length))) {
          
          # define chunk for the target of mapping operation
          .l_current <- purrr::map(.x = .l, .f = ~.x[parallel::splitIndices(nx = length(.l[[1]]), ncl = .no_chunks)[[..i]]])
          
          # define function to be fun
          function_to_run <- function(.l_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
            
            lapply(vector_global_packages, library, character.only = TRUE)
            
            for (i in vector_global_variables) {
              assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
            }
            
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
            
            return(NULL)
            
          }
          
          # assign background worker
          assign(
            x = paste("chunk_", ..i, sep = ""), 
            value = callr::r_bg(
              cmdargs	= c(.globals_save_dir, ..i),
              args = list(".l_current" = .l_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i),
              func = function_to_run
            )
          )
          
          list_workers <- purrr::splice(list_workers, get(paste("chunk_", ..i, sep = "")))
          names(list_workers)[length(list_workers)] <- paste("chunk_", ..i, sep = "")
          
        }
        
        if (.debug == TRUE) {
          global_list_workers <<- list_workers
        }
        
      }
      
      # SPLICER
      
      # initialise list_result if it hasnt already been created
      if (ls(pattern = "^list_result$") %>% length == 0) {
        list_result <- list()
        
        vector_current_chunks_spliced <- numeric()
      }
      
      if (length(vector_current_chunks_spliced) < map_length) {
        
        vector_completed_chunks <- names(vector_logical_indices_workers_completed_reported) %>% gsub(pattern = "chunk_", replacement = "") %>% type.convert(as.is = TRUE)
        
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
        
        Sys.sleep(1)
        
      }
      
      # END SPLICER ###
      
      # an actually useful progress bar although rudimentary
      if (ls(pattern = "^list_result$") %>% length == 0) {
        cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%)", sep = ""))
      } else {
        cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%); Splice progress: ", length(vector_current_chunks_spliced), "/", map_length, sep = ""))
      }
      
      # deal with completion flag
      flag_completion <- (length(vector_current_chunks_spliced) == map_length) & (length(which(vector_logical_indices_workers_completed_reported)) == map_length)
      
    }
    
    vector_exit_statuses <- unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status()))
    
    Sys.sleep(1)
    
  }
  
  if (.keep_intermediate_files == FALSE) {
    # unlink(list.files(path = .intermediate_files_dir, pattern = paste(.job_name, "_chunk_.*.rds", sep = ""), full.names = TRUE ), recursive = TRUE)
    # unlink(list.files(path = .temp_dir, pattern = paste(.job_name, ".rdata", sep = ""), full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .globals_save_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .intermediate_files_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
  }
  
  # if any error, stop all
  if (any(vector_exit_statuses != 0)) {
    purrr::map(.x = list_workers, .f = ~.x$signal(9))
    stop(print(paste("Exit status failure received on chunks:", paste(names(list_workers)[which(vector_exit_statuses != 0)], collapse = ", "))))
  }
  
  rm(list = ls(pattern = "chunk_"))
  rm(list_workers)
  
  cat("\n")
  
  # unchunkify
  list_result <- purrr::flatten(list_result)
  
  return(list_result)
  
  
}