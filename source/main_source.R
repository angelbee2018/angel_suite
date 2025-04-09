# This document contains all the amazing scripts necessary to do regular bioinformatics manipulations.
# PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME 
# For the citation information CFF file, please visit: https://github.com/angel-bee2018/angel_suite/blob/master/CITATION.cff
# PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME PLEASE CITE ME 
## no warranty or guarantee is provided for any failures or losses associated with any use of content in this document. I am just sharing this because somemone else might need something like this, and I've done the hard work already to make it a reality.

# INSIDE-OUT DYNGET
fn_inside_out_dyget <- function(.uncollected_variables, .debug = FALSE, .frame_mark) {
  
  # recursively traverse the frame stack from the bottom up until we collect everything
  temp_current_frame <- sys.nframe()
  
  while (length(.uncollected_variables) > 0 & temp_current_frame > -1) {
    
    for (j in .uncollected_variables) {
      
      assign(x = j, value = dynGet(x = j, ifnotfound = NULL, minframe = temp_current_frame), envir = sys.frame(.frame_mark))
      
      if (.debug == TRUE) {
        message("temp_current_frame: ", temp_current_frame)
        message(j)
        print(dynGet(x = j, ifnotfound = NULL, minframe = .frame_mark))
      }
      
      if (is.null(dynGet(x = j, ifnotfound = NULL, minframe = .frame_mark)) == TRUE) {
        rm(list = j, envir = sys.frame(.frame_mark))
      }
      
    }
    
    .uncollected_variables <- setdiff(.uncollected_variables, ls())
    
    temp_current_frame <- temp_current_frame - 1
    
  }
  
  return(.uncollected_variables)
  
}

# RECURSIVE FUNCTION SCOPING
fn_recursive_function_autoscoper <- function(.vector_functions_for_inspection, .vector_excluded_variables, .debug = FALSE, .frame_mark) {
  
  temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- .vector_functions_for_inspection
  
  temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- character()
  
  # recursively inspect functions until there are no more functions within the functions
  # nested packages -> added to global vector of package names
  # nested functions -> copied into environment using `dynGet` + uncollected function names are added to the uncollected list + added to vector_global_variables + flagged for further internal inspection
  while (length(temp_vector_names_of_newly_discovered_nested_functions_for_inspection) > 0) {
    
    for (temp_function_to_inspect in purrr::map(.x = temp_vector_names_of_newly_discovered_nested_functions_for_inspection, .f = ~dynGet(x = .x, ifnotfound = NULL, minframe = .frame_mark))) {
      
      vector_global_packages <- c(vector_global_packages, globals::packagesOf(globals::globalsOf(temp_function_to_inspect, mustExist = FALSE)) %>% setdiff(., c("base", "rlang"))) %>% unique
      
      temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- c(temp_vector_names_of_newly_discovered_nested_functions_for_inspection0, setdiff(globals::findGlobals(temp_function_to_inspect), .vector_excluded_variables))
      
    }
    
    # flag for further internal inspection
    temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- temp_vector_names_of_newly_discovered_nested_functions_for_inspection0
    
    temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- fn_inside_out_dyget(.uncollected_variables = temp_vector_names_of_newly_discovered_nested_functions_for_inspection0, .frame_mark = .frame_mark)
    
    # uncollected function names are added to the global uncollected list
    .vector_excluded_variables <- c(.vector_excluded_variables, temp_vector_names_of_newly_discovered_nested_functions_for_inspection0) %>% unique
    
    # flag for further internal inspection (update to have only collected and loaded functions)
    temp_vector_names_of_newly_discovered_nested_functions_for_inspection <- intersect(temp_vector_names_of_newly_discovered_nested_functions_for_inspection, lsf.str() %>% as.character)
    
    # collected functions added to vector_global_variables
    vector_global_variables <- c(vector_global_variables, temp_vector_names_of_newly_discovered_nested_functions_for_inspection) %>% unique
    
    temp_vector_names_of_newly_discovered_nested_functions_for_inspection0 <- character()
    
  }
  
  return(.vector_excluded_variables)
  
  
}

# SCOPER
fn_variable_autoscoper <- function(.f, .user_global_objects, .debug = FALSE, .frame_mark, .globals_save_dir) {
  
  vector_global_variables <- c(globals::findGlobals(.f), .user_global_objects)
  
  if (.debug == TRUE) {
    print(vector_global_variables)
  }
  
  temp_global_variables_uncollected <- setdiff(vector_global_variables, ls())
  
  temp_global_variables_uncollected <- fn_inside_out_dyget(.uncollected_variables = temp_global_variables_uncollected, .frame_mark = .frame_mark)
  
  # deal with the special case of nested function definitions
  vector_all_function_names_in_environment <- lsf.str() %>% as.character
  
  temp_global_variables_uncollected <- fn_recursive_function_autoscoper(.vector_functions_for_inspection = vector_all_function_names_in_environment, .vector_excluded_variables = temp_global_variables_uncollected, .frame_mark = .frame_mark)
  
  return(
    list(
      "vector_global_variables" = vector_global_variables,
      "temp_global_variables_uncollected" = temp_global_variables_uncollected
    )
  )
  
}

# simple round robin callr
## A wrapper for `callr:r_bg` that implements parallel computation for `purrr::pmap`. 
## it functions exactly like `furrr` except intermediate results are written to disk, with no chance of memory leakage.
## how it works: <your list> -> automatically inspected for required variables and packages -> variables saved to a temporary file -> r_bg called -> socket transfer of function arguments from parent R process to child R process -> packages automatically loaded in child -> temporary file read back into child -> the real computation proceeds -> list is spliced back together ON THE FLY in the parent process
## this function works best when a) there are large variables (1GB+) referenced in the function, or b) there is extreme memory leakage when using `furrr`/`future`, `biocparallel` or `parallel`, or c) there are recursive multithreading calls inside multithreading child processes, or d) the parent process is slow because there are alot of variables stored in the global environment, or e) the parent process feels like it needs a fresh start but that's not possible.
## against furrr or native pmap, this function won't speed up for: a) tiny computations that don't even need parallel in the first place, b) very simple function bodies, c) input variables which are small in size, d) computations that are already highly compact and efficient. example is Hadley's: boot_df <- function(x) x[sample(nrow(x), replace = T), ]; rsquared <- function(mod) summary(mod)$r.squared; boot_lm <- function(i) {  rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))}

# .l, .f: purrr::pmap arguments
# .num_workers: number of parallel processes to use. this is a maximum value. will never exceed number of chunks. DEFAULT: 1
# chunks: the .l list is divided evenly into chunks. chunks will be calculated in their own parallel process in a round robin fashion.
# .no_chunks: how many chunks do you want, from 1 to the length of .l. DEFAULT: 1
# .splicing_order: affects the way in which the result chunks from workers are stitched back together into a list. two personalities are permitted using the option : 1. ordered - return results in order (at each tick, check if we have consecutive saved files ready before splicing them into a single list); 2. unordered - return list elements as they arrive. 3. any other option: does not return anything, which is useful if you just want to save results to disk. DEFAULT: "ordered"

# NOTE: YOU NEED TO SPECIFY YOUR OWN NUMBER OF CHUNKS AND WORKERS OTHERWISE IT WILL JUST RUN IN SEQUENTIAL. I DON'T KNOW WHAT YOU WANT.
## TYPICAL USAGE: at the bare minimum, just specify .num_workers and .no_chunks and that's all it needs to run!

# .job_name: the name of the job - underscores are automatically added on both sides filenames. DEFAULT: epoch time `as.numeric(Sys.time())` at the time of execution
# .globals_save_path/.temp_path (D): save location of env variables from parent process as .rdata. DEFAULT: paste(tempdir(), "/", as.numeric(Sys.time()), ".rdata", sep = "")
# .globals_save_compress: compress the .rdata object or no. DEFAULT: TRUE
# .re_export: if set to FALSE, if .globals_save_path already exist, then do not recollect and resave env variables from parent process. if TRUE, env variables will always overwrite existing file on disk
# .globals_mode: set this to control how globals are exported from the parent process. "auto": the script will attempt to automatically grab a necessary variables and packages by performing a similar "black magic" ritual as future/furrr and pray it works (oh God help us). "global": parent process will save list = ls().
# .user_global_objects: only for use with "auto" .globals_mode. A character vector of additional variables you want to be exported.
# .intermediate_files_dir/.temp_dir (D): folder location of all intermediate files for the computation. DEFAULT: tempdir()
# .status_messages_dir: folder location of status messages. DEFAULT: tempdir()
# .progress: purrr parameter. DEFAULT: TRUE
# .keep_intermediate_files: TRUE: keep files written to .intermediate_files_dir. FALSE: remove after done. NOTE: this will alway delete files when this script exits gracefully. use .debug to keep intermediate files no matter what.
# .debug: TRUE: keep and report EVERYTHING. DEFAULT: FALSE

round_robin_pmap_callr <- function(
    .l, .f, 
    .num_workers = 1, .no_chunks = 1, .splicing_order = "ordered", 
    .job_name = NULL, 
    .globals_save_dir = NULL, .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL, 
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE, 
    .status_messages_dir = NULL, .progress = TRUE,
    .debug = FALSE) {
  
  # DEBUG ###
  # .l = list(
  #   "b1" = 1:20
  # )
  # .num_workers = 20
  # .globals_mode = "user"
  # .re_export = TRUE
  # .globals_save_dir = paste(tempdir(), "/tempdata.rdata", sep = "")
  # .intermediate_files_dir = tempdir()
  # .user_global_objects = c()
  # .status_messages_dir = paste(tempdir(), sep = "")
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
      # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.globals_save_dir, i, ".qs", sep = ""))
      }
      
    }
    
  }
  
  # END SCOPING ###
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(purrr::map(.x = .l, .f = ~length(.x))))
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  if (.num_workers > map_length) {
    .num_workers <- map_length
  }
  
  if (map_length == 0) {
    return(.l[[1]])
  }
  
  # preallocate worker list
  list_workers <- list()
  
  # CHUNKING ###
  
  if (is.null(.no_chunks) == TRUE) {
    .no_chunks <- .num_workers
  } else if (is.numeric(.no_chunks) == FALSE) {
    .no_chunks <- .num_workers
  }
  
  if (.no_chunks > map_length) {
    .no_chunks <- map_length
  }
  
  map_length <- .no_chunks
  
  # END CHUNKING ###
  
  message("Commencing computation")
  
  # initial allocation of tasks to maximum number of workers
  for (..i in 1:.num_workers) {
    
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
        print(ls(all.names = TRUE))
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
  
  global_list_workers_initial <<- list_workers
  
  # keep track of iteration progress
  current_map_index <- .num_workers
  
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
      if (number_of_processes_alive < .num_workers & current_map_index < map_length) {
        
        new_map_start <- current_map_index + 1
        new_map_end <- min(c(current_map_index + .num_workers - number_of_processes_alive, map_length))
        
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
              print(ls(all.names = TRUE))
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
    
    Sys.sleep(0.1)
    
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

callr_map <- function(
    .x, .f, 
    .no_workers = 1, .no_chunks = 1, .splicing_order = "ordered", 
    .job_name = NULL, 
    .globals_save_dir = NULL, .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL, 
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE, 
    .status_messages_dir = NULL, .progress = TRUE,
    .debug = FALSE) {
  
  # DEBUG ###
  # x = list(1:20)
  # .no_workers = 20
  # .globals_mode = "user"
  # .re_export = TRUE
  # .globals_save_dir = paste(tempdir(), "/tempdata.rdata", sep = "")
  # .intermediate_files_dir = tempdir()
  # .user_global_objects = c()
  # .status_messages_dir = paste(tempdir(), sep = "")
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
    .job_name <- paste("callr_map_", .epoch_time, sep = "")
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
      # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.globals_save_dir, i, ".qs", sep = ""))
      }
      
    }
    
  }
  
  # END SCOPING ###
  
  # check if all the list elements are of equal length
  map_length <- length(.x)
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  if (.no_workers > map_length) {
    .no_workers <- map_length
  }
  
  if (map_length == 0) {
    return(.x)
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
    .x_current <- .x[parallel::splitIndices(nx = length(.x), ncl = .no_chunks)[[..i]]]
    
    # define function to be fun
    function_to_run <- function(.x_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
      
      lapply(vector_global_packages, library, character.only = TRUE)
      
      for (i in vector_global_variables) {
        assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
      }
      
      # load(file = .globals_save_dir, envir = .GlobalEnv)
      
      if (.debug == TRUE) {
        print("ls before map")
        print(ls(all.names = TRUE))
      }
      
      obj <- purrr::map(
        .x = .x_current,
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
        args = list(".x_current" = .x_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i),
        func = function_to_run
      )
    )
    
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
          .x_current <- .x[parallel::splitIndices(nx = length(.x), ncl = .no_chunks)[[..i]]]
          
          # define function to be fun
          function_to_run <- function(.x_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
            
            lapply(vector_global_packages, library, character.only = TRUE)
            
            for (i in vector_global_variables) {
              assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
            }
            
            # load(file = .globals_save_dir, envir = .GlobalEnv)
            
            if (.debug == TRUE) {
              print("ls before map")
              print(ls(all.names = TRUE))
            }
            
            obj <- purrr::map(
              .x = .x_current,
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
              args = list(".x_current" = .x_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i),
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
    
    Sys.sleep(0.1)
    
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

callr_map2 <- function(
    .x, .y, .f, 
    .no_workers = 1, .no_chunks = 1, .splicing_order = "ordered", 
    .job_name = NULL, 
    .globals_save_dir = NULL, .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL, 
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE, 
    .status_messages_dir = NULL, .progress = TRUE,
    .debug = FALSE) {
  
  # DEBUG ###
  # x = list(1:20)
  # .no_workers = 20
  # .globals_mode = "user"
  # .re_export = TRUE
  # .globals_save_dir = paste(tempdir(), "/tempdata.rdata", sep = "")
  # .intermediate_files_dir = tempdir()
  # .user_global_objects = c()
  # .status_messages_dir = paste(tempdir(), sep = "")
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
    .job_name <- paste("callr_map_", .epoch_time, sep = "")
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
      # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.globals_save_dir, i, ".qs", sep = ""))
      }
      
    }
    
  }
  
  # END SCOPING ###
  
  # check if all the list elements are of equal length
  map_length <- unique(c(length(.x), length(.y)))
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  if (.no_workers > map_length) {
    .no_workers <- map_length
  }
  
  if (map_length == 0) {
    return(NULL)
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
    .x_current <- .x[parallel::splitIndices(nx = length(.x), ncl = .no_chunks)[[..i]]]
    .y_current <- .y[parallel::splitIndices(nx = length(.y), ncl = .no_chunks)[[..i]]]
    
    # define function to be fun
    function_to_run <- function(.x_current, .y_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
      
      lapply(vector_global_packages, library, character.only = TRUE)
      
      for (i in vector_global_variables) {
        assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
      }
      
      # load(file = .globals_save_dir, envir = .GlobalEnv)
      
      if (.debug == TRUE) {
        print("ls before map")
        print(ls(all.names = TRUE))
      }
      
      obj <- purrr::map2(
        .x = .x_current,
        .y = .y_current,
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
        args = list(".x_current" = .x_current, ".y_current" = .y_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i),
        func = function_to_run
      )
    )
    
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
          .x_current <- .x[parallel::splitIndices(nx = length(.x), ncl = .no_chunks)[[..i]]]
          .y_current <- .y[parallel::splitIndices(nx = length(.y), ncl = .no_chunks)[[..i]]]
          
          # define function to be fun
          function_to_run <- function(.x_current, .y_current, .f, .globals_save_dir, .intermediate_files_dir, .job_name, .progress, .debug, vector_global_variables, vector_global_packages, ..i) {
            
            lapply(vector_global_packages, library, character.only = TRUE)
            
            for (i in vector_global_variables) {
              assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
            }
            
            # load(file = .globals_save_dir, envir = .GlobalEnv)
            
            if (.debug == TRUE) {
              print("ls before map")
              print(ls(all.names = TRUE))
            }
            
            obj <- purrr::map2(
              .x = .x_current,
              .y = .y_current,
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
              args = list(".x_current" = .x_current, ".y_current" = .y_current, ".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".progress" = .progress, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages, "..i" = ..i),
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
    
    Sys.sleep(0.1)
    
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

# use callr to make an insulated call
callr_insulator <- function(
    .f,
    .num_workers = 1, .no_chunks = 1, 
    .job_name = NULL,
    .globals_save_dir = NULL, .globals_save_compress = TRUE, .re_export = TRUE, .globals_mode = "auto", .user_global_objects = NULL,
    .intermediate_files_dir = NULL, .keep_intermediate_files = FALSE,
    .status_messages_dir = NULL,
    .debug = FALSE) {
  
  # DEBUG ###
  # .f = combn(x = 1:8, m = 3)
  # .num_workers = 1
  # .no_chunks = 1
  # .job_name = NULL
  # .globals_save_dir = NULL
  # .globals_save_compress = TRUE
  # .re_export = TRUE
  # .globals_mode = "auto"
  # .user_global_objects = NULL
  # .intermediate_files_dir = NULL
  # .keep_intermediate_files = FALSE
  # .status_messages_dir = NULL
  # .debug = FALSE
  ###########

  if (is.call(.f) == FALSE & (is.function(.f) == FALSE)) {
        stop("Expect .f to be of class \"call\" or \"function\"")
        }
  
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
  
  # # SCOPING - FUNCTIONED BUT CURRENTLY IS BROKEN ###
  # 
  # vector_global_packages <- globals::packagesOf(globals::globalsOf(.f, mustExist = FALSE)) %>% setdiff(., c("base", "rlang"))
  # 
  # if (.debug == TRUE) {
  #   print("vector_global_packages at the surface level of the function")
  #   print(vector_global_packages)
  # }
  # 
  # ## mark the current frame number
  # .frame_mark <- sys.nframe()
  # 
  # list_fn_variable_autoscoper <- fn_variable_autoscoper(.f = .f, .user_global_objects = .user_global_objects, .debug = .debug, .frame_mark = .frame_mark, .globals_save_dir = .globals_save_dir)
  # 
  # vector_global_variables <- list_fn_variable_autoscoper$vector_global_variables
  # temp_global_variables_uncollected <- list_fn_variable_autoscoper$temp_global_variables_uncollected
  # 
  # if (length(temp_global_variables_uncollected) > 0) {
  #   warning(paste("Scoping has finished but there remain some uncollected variables: ", paste(temp_global_variables_uncollected, collapse = " , "), "\n"))
  #   
  #   if (.debug == TRUE) {
  #     message("all variables in environment right before writing to disk")
  #     print(ls())
  #   }
  #   
  #   vector_global_variables <- setdiff(vector_global_variables, temp_global_variables_uncollected)
  # }
  # 
  # message("Writing globals to disk")
  # # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
  # 
  # for (i in vector_global_variables) {
  #   qs::qsave(x = get(i, envir = sys.frame(.frame_mark)), file = paste(.globals_save_dir, i, ".qs", sep = ""))
  # }
  # 
  # # END SCOPING - FUNCTIONED ###
  
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
        
        if (.debug == TRUE) {
          print("ls_1")
          print(ls(all.names = TRUE))
        }
        
        temp_global_variables_uncollected <- setdiff(temp_global_variables_uncollected, ls(all.names = TRUE))
        
        if (.debug == TRUE) {
          print("temp_global_variables_uncollected")
          print(temp_global_variables_uncollected)
        }
        
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
      # save(list = vector_global_variables, file = .globals_save_dir, compress = .globals_save_compress)
      
      for (i in vector_global_variables) {
        qs::qsave(x = get(i), file = paste(.globals_save_dir, i, ".qs", sep = ""))
      }
      
    }
    
  }
  
  # END SCOPING ###
  
  message("Commencing insulated computation")
  
  # define function to be fun
  function_to_run <- function(.f, .globals_save_dir, .intermediate_files_dir, .job_name, .debug, vector_global_variables, vector_global_packages) {
    
    lapply(vector_global_packages, library, character.only = TRUE)
    
    for (i in vector_global_variables) {
      assign(x = i, value = qs::qread(file = paste(.globals_save_dir, i, ".qs", sep = "")), envir = .GlobalEnv)
    }
    
    # load(file = .globals_save_dir, envir = .GlobalEnv)
    
    if (.debug == TRUE) {
      print("ls before pmap")
      print(ls(all.names = TRUE))
    }

    if (is.call(.f) == TRUE) {
      obj <- eval(.f)      
      } else if (is.function(.f) == TRUE) {
      obj <- .f()
      }
    
    # saveRDS(object = obj, file = paste(.intermediate_files_dir, "/", .job_name, "_chunk_", ..i, ".rds", sep = ""))
    qs::qsave(x = obj, file = paste(.intermediate_files_dir, "/", .job_name, "_insulator.qs", sep = ""))
    
    return(NULL)
    
  }
  
  # assign background worker
  assign(
    x = "callr_insulator",
    value = callr::r_bg(
      cmdargs	= c(.globals_save_dir, "callr_insulator"),
      args = list(".f" = .f, ".globals_save_dir" = .globals_save_dir, ".intermediate_files_dir" = .intermediate_files_dir, ".job_name" = .job_name, ".debug" = .debug, "vector_global_variables" = vector_global_variables, "vector_global_packages" = vector_global_packages),
      func = function_to_run
    )
  )
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  flag_completion <- FALSE

  while(flag_completion == FALSE) {
    
    # check on process status and write progress file
    stdout_lines <- callr_insulator$read_output_lines()
    stderr_lines <- callr_insulator$read_error_lines()
    
    if (length(stdout_lines) != 0) {
      write(x = paste(Sys.time()), file = paste(.status_messages_dir, .job_name, "_", "callr_insulator", "_stdout.txt", sep = ""), append = TRUE)
      write.table(x = stdout_lines, file = paste(.status_messages_dir, .job_name, "_", "callr_insulator", "_stdout.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE) %>% suppressWarnings()
    }
    
    if (length(stderr_lines) != 0) {
      write(x = paste(Sys.time()), file = paste(.status_messages_dir, .job_name, "_", "callr_insulator", "_stderr.txt", sep = ""), append = TRUE)
      write.table(x = stderr_lines, file = paste(.status_messages_dir, .job_name, "_", "callr_insulator", "_stderr.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE) %>% suppressWarnings()
    }
    
    # if any error, stop all
    insulator_exit_status <- callr_insulator$get_exit_status()
    
    if (!is.null(insulator_exit_status)) {
      if (insulator_exit_status != 0) {
        callr_insulator$signal(9)
        stop(print("Exit status failure received in insulator"))
      }
    }
    
    # check for errors that don't show up in process exit status
    logical_any_error <- grepl(x = stderr_lines, pattern = "^Error in", ignore.case = FALSE)
    
    if (length(logical_any_error) > 0) {
      if (any(logical_any_error) == TRUE) {
        callr_insulator$signal(9)
        stop(print("Stop error received"))
      }
    }
    
    # completion flag
    flag_isalive <- callr_insulator$is_alive()
    flag_exit_status_is_zero <- if (is.null(insulator_exit_status)) {FALSE} else {insulator_exit_status == 0}
    flag_completion <- (flag_isalive == FALSE) & (flag_exit_status_is_zero == TRUE)
    
    if (flag_completion == TRUE) {
      break()
    } else if (flag_completion == FALSE) {
      Sys.sleep(0.1)
    }
    
  }
  
  result <- qs::qread(file = paste(.intermediate_files_dir, "/", .job_name, "_insulator.qs", sep = ""))
  
  if (.keep_intermediate_files == FALSE) {
    # unlink(list.files(path = .intermediate_files_dir, pattern = paste(.job_name, "_chunk_.*.rds", sep = ""), full.names = TRUE ), recursive = TRUE)
    # unlink(list.files(path = .temp_dir, pattern = paste(.job_name, ".rdata", sep = ""), full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .temp_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
    unlink(list.files(path = .intermediate_files_dir, pattern = ".*.qs", full.names = TRUE ), recursive = TRUE)
  }
  
  rm(callr_insulator)
  
  cat("\n")
  
  return(result)
  
}

# test <- round_robin_pmap_callr(
#   .l = list(
#     "b1" = a1 %>% purrr::array_tree() %>% head,
#     "b2" = 1:length(a1 %>% purrr::array_tree() %>% head)
#   ),
#   .num_workers = 4,
#   .globals_mode = "user",
#   .re_export = TRUE,
#   .globals_save_path = paste(tempdir, "/tempdata.rdata", sep = ""),
#   .intermediate_files_dir = tempdir,
#   .user_global_objects = c(),
#   .status_messages_dir = paste(tempdir, sep = ""),
#   .job_name = "test",
#   .splicing_order = "unordered",
#   .f = function(b1, b2) {Sys.sleep(0); return(list(LETTERS[b2]))})

# function to aggregate a bunch of individual numbers into discrete intervals
numbers_to_intervals <- function(vector_input) {
  
  vector_input_sorted <- sort(unlist(vector_input))
  
  # convert the individual nt positions to a set of ranges
  ## achieve this by comparing n to n + 1
  ## for rows where diff > 1, the n represents the end of an exon. n.plus.1 represents the start of the exon right after the gap.
  if (length(vector_input_sorted) == 1) {
    
    vec_starts <- vector_input_sorted
    vec_ends <- vector_input_sorted
    vec_widths <- 1
    
  } else {
    
    tibble_n_n.plus.1 <- tibble::tibble("n" = head(x = vector_input_sorted, n = length(vector_input_sorted) - 1),
                                        "n.plus.1" = vector_input_sorted[2:length(vector_input_sorted)])
    
    tibble_n_n.plus.1 <- tibble::add_column(.data = tibble_n_n.plus.1, "difference" = tibble_n_n.plus.1$n.plus.1 - tibble_n_n.plus.1$n)
    
    ## if there are no gaps, then just take the range as the start:end
    if (length(unique(tibble_n_n.plus.1$difference)) == 1) {
      
      vec_starts <- min(vector_input_sorted)
      vec_ends <- max(vector_input_sorted)
      vec_widths <- vec_ends - vec_starts + 1
      
    } else if (length(unique(tibble_n_n.plus.1$difference)) > 1) {
      
      # vec_blockCount <- tibble_n_n.plus.1$difference %>% unique %>% length
      vec_starts <- c(vector_input_sorted[1], 
                      unlist(tibble_n_n.plus.1[tibble_n_n.plus.1$difference > 1, "n.plus.1"], use.names = FALSE))
      vec_ends <- c(unlist(tibble_n_n.plus.1[tibble_n_n.plus.1$difference > 1, "n"], use.names = FALSE), 
                    vector_input_sorted[length(vector_input_sorted)])
      vec_widths <- vec_ends - vec_starts + 1
      
    }
    
  }
  
  tibble_intervals <- tibble::tibble(
    "start" = vec_starts,
    "end" = vec_ends,
    "width" = vec_widths 
  )
  
  return(tibble_intervals)
  
}

# function to do Benjamini-Hochberg FDR correction for Gene Ontology output tables from GOHyperGall

## takes the output of the GOHyperGall GO enrichment table and over-writes the Padj column with benjamini-corrected values from the Phyper column
## spits out the resulting table with modified single column.

## NOTE: BENJAMINI PVALUES ABOVE 0.05 ARE RENAMED NA

GOHyperGAll_benjamini_correction <- function(raw_GOHyperGAll_table)
  
{
  
  benjamini_GOHyperGAll_table <- raw_GOHyperGAll_table
  
  benjamini_GOHyperGAll_table <- benjamini_GOHyperGAll_table[benjamini_GOHyperGAll_table$SampleMatch > 1, ]
  
  if (benjamini_GOHyperGAll_table %>% nrow == 0) {
    benjamini_GOHyperGAll_table <- tibble()
  } else {
    benjamini_GOHyperGAll_table[, "Padj"] <- p.adjust(p = benjamini_GOHyperGAll_table[, "Phyper"], method = "BH", n = length(benjamini_GOHyperGAll_table[, "Phyper"]))
    
    # benjamini_GOHyperGAll_table[benjamini_GOHyperGAll_table$Padj > 0.05, "Padj"] <- NA
  }
  
  return(benjamini_GOHyperGAll_table)
  
}

# the equivalent for bc3net::enrichment() output table

bc3net_benjamini_correction <- function(raw_bc3net_table)
  
{
  
  benjamini_bc3net_table <- raw_bc3net_table
  
  benjamini_bc3net_table_processed <- benjamini_bc3net_table[benjamini_bc3net_table$genes > 1, ]
  
  if (nrow(benjamini_bc3net_table_processed) != 0) {
    
    benjamini_bc3net_table_processed[, "padj"] <- p.adjust(p = benjamini_bc3net_table_processed[, "pval"], method = "BH", n = length(benjamini_bc3net_table_processed[, "pval"]))
    
    # benjamini_bc3net_table_processed[benjamini_bc3net_table_processed$padj > 0.05, "padj"] <- NA
    
  } else {
    
    benjamini_bc3net_table_processed <- benjamini_bc3net_table
    
    benjamini_bc3net_table_processed[, "padj"] <- NA
    
  }
  
  return(benjamini_bc3net_table_processed)
  
}

# DEPRECATED CRYPTIC FUNCTION FOR BC3NET CATALOG GENERATION
## bc3net::enrichment() does not show captured genes for each family enriched, so we have to add it in. but in doing so, i want to avoid a purrr within a purrr

## this function selects genes from the background in each family which are ONLY input genes.

# filtering_genehits_from_background_catalogue <- function(catalogue, genehit_vector){
#   
#   filtered_catalogue <- purrr::map(.x = catalogue, .f = ~intersect(.x, genehit_vector))
#   
#   return(filtered_catalogue)
#   
# }

# 

## FUNCTION TO EXTRACT REFERENCE EXONS WHICH OVERLAP EXACTLY WITH QUERY EXONS
## NOTE: to be used with purrr
## input: spliceregion_list: a list containing details of ONE junction: $chr, $diff_exon_start, $diff_exon_end
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_matching.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type = "exon") {
  
  # DEBUG ###################
  # index <- 1
  # spliceregion_list <- wide_tibble_of_all_unique_VSR_and_exon_coords_array.tree_not_IR[[index]]
  # # tibble_gtf_table <- tibble_ref_gtf
  # tibble_gtf_table <- tibble_recon_gtf
  # stranded = FALSE
  ###########################
  
  # print(query_chr, "\n")
  # print(query_start, "\n")
  # print(query_end, "\n")
  # print(query_strand, "\n")
  # print(tibble_gtf_table, "\n")
  
  # print(paste("now processing junction number", index))
  
  global_query_chr <<- query_chr
  global_query_start <<- query_start
  global_query_end <<- query_end
  global_query_strand <<- query_strand
  
  if (is.na(query_strand) | is.null(query_strand)) {
    query_strand <- "*"
  }
  
  if (!(query_strand == "+" | query_strand == "-")) {
    
    # cat("a\n")
    
    # +/- 1 nt tolerance
    tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
      .[which(.$start > ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance) & .$end < ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance)), ] %>% 
      .[which((.$start < ((query_start %>% as.numeric) + 1 + left_query_shift + left_tolerance) & .$end > ((query_end %>% as.numeric) - 1 + right_query_shift - right_tolerance))), ] %>% 
      .[which(.$type %in% return_type), ]
    
  } else if (query_strand == "+" | query_strand == "-") {
    
    # cat("b\n")
    
    # +/- 1 nt tolerance
    tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                           tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
      .[which(.$start > ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance) & .$end < ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance)), ] %>% 
      .[which((.$start < ((query_start %>% as.numeric) + 1 + left_query_shift + left_tolerance) & .$end > ((query_end %>% as.numeric) - 1 + right_query_shift - right_tolerance))), ] %>% 
      .[which(.$type %in% return_type), ]
    
  }
  
  return(tibble_gtf_subset_matching_exons)
  
}

## END extract_matching.exons() ###

## FUNCTION TO EXTRACT TRANSCRIPTS WITH JUNCTION-FLANKING EXONS.
## NOTE: to be used with purrr
## details of ONE junction: $chr, $start, $end, $strand
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_junction.flanking.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, match_consecutive = TRUE, match_inside_same_transcript = TRUE, return_type = "exon") {
  
  # DEBUG ###################
  
  # query_chr = query_chr
  # query_start = a1$event_region_start %>% type.convert
  # query_end = a1$event_region_end %>% type.convert
  # query_strand = "*"
  # tibble_gtf_table = tibble_ref_gtf
  # tolerance_left = 0
  # tolerance_right = 0
  # tolerance_inside = 0
  # tolerance_outside = 0
  # match_consecutive = FALSE
  # return_type = "exon"
  
  ###########################
  
  # print(paste("now processing junction number", index))
  
  # global_query_chr <<- query_chr
  # global_query_start <<- query_start
  # global_query_end <<- query_end
  # global_query_strand <<- query_strand
  
  if (is.na(query_strand) | is.null(query_strand)) {
    query_strand <- "*"
  }
  
  if (query_strand == "." | query_strand == 0 | query_strand == "*") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance) & .$end >= ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance), ] %>% .[!(.$start <= ((query_end %>% as.numeric) + right_query_shift - right_tolerance) & .$end >= ((query_start %>% as.numeric) + left_query_shift + left_tolerance)), ] %>% .[.$type %in% return_type, ]
    
  } else if (query_strand == "+" | query_strand == "-") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$strand == query_strand %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + right_query_shift + right_tolerance) & .$end >= ((query_start %>% as.numeric) - 1 + left_query_shift - left_tolerance), ] %>% .[!(.$start <= ((query_end %>% as.numeric) + right_query_shift - right_tolerance) & .$end >= ((query_start %>% as.numeric) + left_query_shift + left_tolerance)), ] %>% .[.$type %in% return_type, ]
    
  }
  
  list_of_junction_associated_transcripts <- tibble_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks a junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  
  if (match_inside_same_transcript == FALSE) {
    
    list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts)
    
  } else if (match_inside_same_transcript == TRUE & match_consecutive == TRUE) {
    
    list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
    
  } else if (match_inside_same_transcript == TRUE & match_consecutive == FALSE) {
    
    list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2)
    
  }
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

## END extract_junction.flanking.exons() ###

## FUNCTION TO EXTRACT REFERENCE EXONS WHICH OVERLAP EXACTLY WITH QUERY EXONS
## NOTE: to be used with purrr
## input: spliceregion_list: a list containing details of ONE junction: $chr, $diff_exon_start, $diff_exon_end
## tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
## index: loop progress marker to be used with imap

extract_overlapping_features <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type = NULL, complete_overlap = FALSE) {
  
  # DEBUG ###################
  # index <- 1
  # spliceregion_list <- wide_tibble_of_all_unique_VSR_and_exon_coords_array.tree_not_IR[[index]]
  # # tibble_gtf_table <- tibble_ref_gtf
  # tibble_gtf_table <- tibble_recon_gtf
  # stranded = FALSE
  ###########################
  
  # cat(query_chr, "\n")
  # cat(query_start, "\n")
  # cat(query_end, "\n")
  # cat(query_strand, "\n")
  # print(tibble_gtf_table, "\n")
  # cat(left_query_shift, "\n")
  # cat(right_query_shift, "\n")
  # cat(left_tolerance, "\n")
  # cat(right_tolerance, "\n")
  
  # print(paste("now processing junction number", index))
  
  global_query_chr <<- query_chr
  global_query_start <<- query_start
  global_query_end <<- query_end
  global_query_strand <<- query_strand
  
  if (is.na(query_strand) | is.null(query_strand)) {
    query_strand <- "*"
  }
  
  if (complete_overlap == FALSE) {
    
    if (!(query_strand == "+" | query_strand == "-")) {
      
      # +/- 1 nt tolerance
      tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
        .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ]
      
    } else if (query_strand == "+" | query_strand == "-") {
      
      # +/- 1 nt tolerance
      tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                             tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
        .[which(.$end >= ((query_start %>% as.numeric) + left_query_shift - left_tolerance) & .$start <= ((query_end %>% as.numeric) + right_query_shift + right_tolerance)), ]
      
    }
    
  } else {
    
    if (!(query_strand == "+" | query_strand == "-")) {
      
      # +/- 1 nt tolerance
      tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
        .[which(.$start <= ((query_start %>% as.numeric) + left_query_shift + left_tolerance) & .$end >= ((query_end %>% as.numeric) + right_query_shift - right_tolerance)), ]
      
    } else if (query_strand == "+" | query_strand == "-") {
      
      # +/- 1 nt tolerance
      tibble_gtf_subset_matching_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                             tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
        .[which(.$start <= ((query_start %>% as.numeric) + left_query_shift + left_tolerance) & .$end >= ((query_end %>% as.numeric) + right_query_shift - right_tolerance)), ]
      
    }
    
  }
  
  if (is.null(return_type) == FALSE) {
    return(tibble_gtf_subset_matching_exons[tibble_gtf_subset_matching_exons$type %in% return_type, ])
  } else if (is.null(return_type) == TRUE) {
    return(tibble_gtf_subset_matching_exons)
  }
  
  return(tibble_gtf_subset_matching_exons)
  
}

## END extract_overlapping_features() ###

## FUNCTIONS TO MAGNETISE A USER VALUE TO REFERENCE START COORDS
### ref_start/end_shift: use +/-1 to get intronic starts/ends
### match to reference starts
magnetise_genome_position_to_ref_starts <- function(query_chr, query_coord, query_strand, tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_start_shift = 0, return_type) {
  
  # DEBUG ###################
  # query_chr <- query_chr
  # query_coord <- 153411582
  # query_strand <- query_strand
  # query_shift <- right_query_shift
  # query_tolerance <- right_tolerance
  # ref_start_shift <- -1
  # return_type <- "exon"
  ###########################
  
  global_magnetised_query_chr <<- query_chr
  global_magnetised_query_coord <<- query_coord
  global_magnetised_query_strand <<- query_strand
  
  if (is.na(query_strand) | is.null(query_strand)) {
    query_strand <- "*"
  }
  
  if (!(query_strand == "+" | query_strand == "-")) {
    
    # +/- 1 nt tolerance
    tibble_ref_starts_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
      .[which((.$start + ref_start_shift) >= (query_coord + query_shift - query_tolerance) & (.$start + ref_start_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
      .[which(.$type %in% return_type), ]
  } else if (query_strand == "+" | query_strand == "-") {
    
    # +/- 1 nt tolerance
    tibble_ref_starts_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws & tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
      .[which((.$start + ref_start_shift) >= (query_coord + query_shift - query_tolerance) & (.$start + ref_start_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
      .[which(.$type %in% return_type), ]
    
  }
  
  magnetised_coord <- (tibble_ref_starts_matched_to_query_coord$start + ref_start_shift) %>%
    .[which(abs(tibble_ref_starts_matched_to_query_coord$start - query_coord) == min(abs(tibble_ref_starts_matched_to_query_coord$start - query_coord)))] %>% 
    .[1]
  
  return(list(
    "tibble_ref_starts_matched_to_query_coord" = tibble_ref_starts_matched_to_query_coord,
    "magnetised_coord" = magnetised_coord
  ) )
  
}

## END magnetise_genome_position_to_ref_starts() ###

### match to reference ends
magnetise_genome_position_to_ref_end <- function(query_chr, query_coord, query_strand, tibble_gtf_table, query_shift = 0, query_tolerance = 0, ref_end_shift = 0, return_type) {
  
  # DEBUG ###################
  # query_chr <- query_chr
  # query_coord <- query_VSR_end
  # query_strand <- query_strand
  # query_shift <- right_query_shift
  # query_tolerance <- right_tolerance
  # ref_end_shift <- -1
  ###########################
  
  if (!(query_strand == "+" | query_strand == "-")) {
    
    # +/- 1 nt tolerance
    tibble_ref_ends_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
      .[which((.$end + ref_end_shift) >= (query_coord + query_shift - query_tolerance) & (.$end + ref_end_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
      .[which(.$type %in% return_type), ]
  } else if (query_strand == "+" | query_strand == "-") {
    
    # +/- 1 nt tolerance
    tibble_ref_ends_matched_to_query_coord <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws & tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
      .[which((.$end + ref_end_shift) >= (query_coord + query_shift - query_tolerance) & (.$end + ref_end_shift) <= (query_coord + query_shift + query_tolerance)), ] %>% 
      .[which(.$type %in% return_type), ]
    
  }
  
  magnetised_coord <- (tibble_ref_ends_matched_to_query_coord$end + ref_end_shift) %>%
    .[which(abs(tibble_ref_ends_matched_to_query_coord$end - query_coord) == min(abs(tibble_ref_ends_matched_to_query_coord$end - query_coord)))] %>% 
    .[1]
  
  return(list(
    "tibble_ref_ends_matched_to_query_coord" = tibble_ref_ends_matched_to_query_coord,
    "magnetised_coord" = magnetised_coord
  ) )
  
}

## END magnetise_genome_position_to_ref_starts() ###

# FUNCTION TO EXTRACT COMMON STRING FROM TWO INPUT STRINGS IN ONE STEP, A SIMPLIFICATION OF QUALV

extract_common_string <- function(input_string_a, input_string_b) {
  
  # debug ###
  # 
  # input_string_a <- "ud_absolute.psi_1"
  # input_string_b <- "ud_absolute.psi_2"
  
  ###########
  
  vector_of_letters_a <- strsplit(input_string_a, "") %>% unlist
  vector_of_letters_b <- strsplit(input_string_b, "") %>% unlist 
  
  raw_LCS <- qualV::LCS(vector_of_letters_a, vector_of_letters_b)
  vector_common_string <- raw_LCS$LCS %>% paste(collapse = "")
  
  return(vector_common_string)
  
}

# END extract_common_string()

# FUNCTION TO TAKE THE AVERAGE VALUE OF A MATRIX OF TIMEPOINTS WITH THREE REPLICATES EACH (3 column compartments at a time)
## replicates of the same timepoint must be all grouped together.
calculate_average_values_from_replicate_columns <- function(input_matrix, number_of_replicates, append_average_to_column_name = TRUE) {
  
  # DEBUG ######
  
  # input_matrix <- wide_tibble_of_psisigma_results_allcomparisons_final_ud.only[, col_indices_observations]
  # number_of_replicates <- 3
  
  ##############
  
  # sanity check - if the number of columns is not an integer multiple of the number of replicates, there's something wrong
  if (ncol(input_matrix) %% number_of_replicates != 0) {
    
    stop("the number of columns in the matrix is not an integer multiple of the number of replicates specified. please check the matrix and try again.")
    
  }
  
  # get the row numbers to subset
  ## start of each compartment
  a <- seq(1, (ncol(input_matrix) - number_of_replicates + 1), number_of_replicates)
  ## end of each compartment
  b <- seq(number_of_replicates, ncol(input_matrix), number_of_replicates)
  # create list of compartments
  c <- purrr::map2(.x = a, .y = b, .f = ~.x:.y)
  # map each column into each compartment
  d <- purrr::map(.x = c, .f = ~input_matrix[, .x])
  # apply the average
  e <- purrr::map(.x = d, .f = ~apply(X = .x, MARGIN = 1, FUN = function(X){mean(X, na.rm = TRUE)}))
  # retrieve the column names of each compartment
  column_names <- purrr::map(.x = d, .f = ~colnames(.x) %>% purrr::reduce(extract_common_string))
  
  if (append_average_to_column_name == TRUE) {
    
    column_names <- paste(column_names, "average", sep = "")
    
  }
  
  # reframe into tibble, return
  f <- e %>% set_names(column_names) %>% as_tibble
  
}

# END calculate_average_values_from_replicate_columns()

# target_colnames: a vector of names of the columns to split at the same time
split_delimited_columns_in_table <- function(input_table, target_colnames, split, columns_to_deduplicate = NULL, .parallel = "none", .no_workers = 1, .no_chunks = 1, .progress = TRUE) {
  
  # DEBUG ###
  # input_table <- tibble("a" = c("345,345", "asdf,5t345", "32454rtg,54", "3245345,rr"), "b" = c("345,345", "asdf,5t345", "32454rtg,54", "3245345,rr"), "c" = "haha")
  # target_colnames <- c("a", "b")
  # split = "\\,"
  # columns_to_deduplicate <- "c"
  ###########
  
  ## check for length equivalence
  ## if multiple columns were selected for splitting and any row has unequal lengths after strsplit then we cannot continue.
  ## using `apply` in this way will return a list of lists, with L1 being columns and L2 being rows of each column
  list_split_columns <- apply(X = input_table[, target_colnames], MARGIN = 2, FUN = function(x) {return(strsplit(x, split = split))} )
  
  if (.parallel == "none") {
    list_split_lengths <- pmap(
      .l = list_split_columns, 
      .f = function(...) {
        return(unique(unlist(purrr::map(.x = list(...), .f = ~length(.x)))))
      },
      .progress = .progress
    )
  } else if (.parallel == "callr") {
    list_split_lengths <- round_robin_pmap_callr(
      .l = list_split_columns, 
      .num_workers = .no_workers, .no_chunks = .no_chunks,
      .f = function(...) {
        return(unique(unlist(purrr::map(.x = list(...), .f = ~length(.x)))))
      },
      .progress = .progress
    )
  } else if (.parallel == "mc") {
    list_split_lengths <- mc_pmap(
      .l = list_split_columns, 
      .no_workers = .no_workers, .no_chunks = .no_chunks,
      .f = function(...) {
        return(unique(unlist(purrr::map(.x = list(...), .f = ~length(.x)))))
      },
      .progress = .progress
    )
  }
  
  if (unique(unlist(purrr::map(.x = list_split_lengths, .f = ~length(.x))))[1] != 1 | length(unique(unlist(purrr::map(.x = list_split_lengths, .f = ~length(.x))))) != 1) {
    stop("Multicolumn splits are uneven. Cannot proceed further.")
  }
  
  # repeat table according to the split lengths
  row_indices_of_table_repeated_by_split <- rep(x = 1:nrow(input_table), times = unlist(list_split_lengths))
  
  input_table_repeated_by_split <- input_table[row_indices_of_table_repeated_by_split, ]
  
  # replace target column with split values
  input_table_repeated_by_split[, target_colnames] <- list_split_columns %>% purrr::map(.f = ~.x %>% unlist) %>% tibble::as_tibble()
  
  split_table <- input_table_repeated_by_split
  
  # if specified, append an index to a particular column
  if (is.null(columns_to_deduplicate) == FALSE) {
    
    # get the duplicated row indices where split lengths > 1
    indices_of_duplicates <- which(unlist(list_split_lengths) > 1)
    
    # get the repetition number where split lengths > 1
    repetition_numbers_of_duplicates <- unlist(list_split_lengths)[which(unlist(list_split_lengths) > 1)]
    
    # list-ify the columns to be appended
    list_deduplicated_columns <- input_table[, columns_to_deduplicate] %>% array_tree(margin = 2) %>% purrr::map(~array_tree(.x))
    
    # map over each column, split the target element and add _[0-9]+
    if (.parallel == "none") {
      list_deduplicated_columns_split <- purrr::map(
        .x = list_deduplicated_columns, 
        .f = function(a1) {
          
          # map a subset each of the L2 (elements of a column)
          a1[indices_of_duplicates] <- purrr::map2(.x = a1[indices_of_duplicates], .y = repetition_numbers_of_duplicates, 
                                                   .f = ~rep(.x, times = .y) %>% unlist %>% paste(., 1:.y, sep = "_"))
          
          return(a1 %>% unlist)
          
        },
        .progress = .progress )
    } else if (.parallel == "callr") {
      list_split_lengths <- round_robin_pmap_callr(
        .l = list(
          "a1" = list_deduplicated_columns
        ), 
        .num_workers = .no_workers, .no_chunks = .no_chunks,
        .f = function(a1) {
          
          # map a subset each of the L2 (elements of a column)
          a1[indices_of_duplicates] <- purrr::map2(.x = a1[indices_of_duplicates], .y = repetition_numbers_of_duplicates, 
                                                   .f = ~rep(.x, times = .y) %>% unlist %>% paste(., 1:.y, sep = "_"))
          
          return(a1 %>% unlist)
          
        },
        .progress = .progress
      )
    } else if (.parallel == "mc") {
      list_split_lengths <- mc_pmap(
        .l = list(
          "a1" = list_deduplicated_columns
        ), 
        .no_workers = .no_workers, .no_chunks = .no_chunks,
        .f = function(a1) {
          
          # map a subset each of the L2 (elements of a column)
          a1[indices_of_duplicates] <- purrr::map2(.x = a1[indices_of_duplicates], .y = repetition_numbers_of_duplicates, 
                                                   .f = ~rep(.x, times = .y) %>% unlist %>% paste(., 1:.y, sep = "_"))
          
          return(a1 %>% unlist)
          
        },
        .progress = .progress
      )
    }
    
    # tibblise
    tibble_deduplicated_columns_split <- list_deduplicated_columns_split %>% as_tibble
    
    # add back every row onto the split table
    for (dedupe_colname in columns_to_deduplicate) {
      
      split_table[, dedupe_colname] <- tibble_deduplicated_columns_split[, dedupe_colname]
      
    }
    
  }
  
  return(split_table %>% type_convert)
  
}

# END split_delimited_column_in_table()

# create new, faster script based on delimiter counts with stringr::str_count
# target_colnames: a vector of names of the columns to split at the same time
# split: a regular expression
split_delimited_columns_in_table2 <- function(input_table, target_colnames, split, columns_to_deduplicate = NULL) {
  
  # DEBUG ###
  # input_table <- tibble("a" = c("345,345", "asdf,5t345", "32454rtg,54", "3245345,rr"), "b" = c("345,345", "asdf,5t345", "32454rtg,54", "3245345,rr"), "c" = "haha")
  # target_colnames <- c("a", "b")
  # split = "\\,"
  # columns_to_deduplicate <- "c"
  ###########
  
  # use stringr::str_count to count the number of desired delimiters in each row
  nummatrix_split_lengths <- apply(X = input_table[, target_colnames], MARGIN = 2, FUN = function(x) {return(stringr::str_count(string = x, pattern = split))} )
  
  # check for length equivalence
  nummatrix_split_lengths_unique <- t(unique(t(nummatrix_split_lengths))) + 1
  number_of_unique_column_split_structures <- ncol(nummatrix_split_lengths_unique)
  
  if (number_of_unique_column_split_structures != 1) {
    stop(paste("Multicolumn splits are uneven. Cannot proceed further. Got number_of_unique_column_split_structures = ", number_of_unique_column_split_structures, sep = ""))
  }
  
  # repeat table according to the split lengths
  vector_row_indices_of_table_repeated_by_split <- rep(x = 1:nrow(input_table), times = nummatrix_split_lengths_unique)
  
  input_table_repeated_by_split <- input_table[vector_row_indices_of_table_repeated_by_split, ]
  
  # replace target column with split values
  list_split_columns <- apply(X = input_table[, target_colnames], MARGIN = 2, FUN = function(x) {return(strsplit(x, split = split))} )
  
  input_table_repeated_by_split[, target_colnames] <- list_split_columns %>% purrr::map(.f = ~.x %>% unlist) %>% tibble::as_tibble()
  
  split_table <- input_table_repeated_by_split
  
  # if specified, append an index to a particular column
  if (is.null(intersect(colnames(input_table), columns_to_deduplicate)) == FALSE) {
    
    vector_repetition_mask <- rep(x = nummatrix_split_lengths_unique, times = nummatrix_split_lengths_unique)
    vector_repetition_mask[vector_repetition_mask == 1] <- 0
    vector_max_repetitions <- unique(vector_repetition_mask) %>% .[. != 0]
    
    for (temp_max_repetition in vector_max_repetitions) {
      vector_repetition_mask[which(vector_repetition_mask == vector_max_repetitions[1])] <- seq(from = 1, to = temp_max_repetition)
    }
    
    vector_repetition_mask <- paste("_", vector_repetition_mask, sep = "")
    
    vector_repetition_mask[vector_repetition_mask == "_0"] <- ""
    
    # add back every row onto the split table
    for (temp_column_to_deduplicate in intersect(colnames(input_table), columns_to_deduplicate)) {
      split_table[, temp_column_to_deduplicate] <- paste(split_table[[temp_column_to_deduplicate]], vector_repetition_mask, sep = "")
    }
    
  }
  
  split_table <- readr::type_convert(split_table)
  
  return(split_table)
  
}

# END split_delimited_columns_in_table2()

# FUNCTION TO PLOT UMAP FOR SAMPLE AND REPLICATE
# tibble_metadata: this should have n columns equal to the number of samples, and m rows equal to number of features (including replicates). basically a collection of vectors of length equal to number of samples.
# all input tables i write will expect that columns are samples and rows are features. any packages that expect the transpose will be handled by the wrapper.
# protected colname: "condition_names", "replicate_names"
plot_UMAP_sample_and_replicate_plotly <- function(table_matrix, condition_order = NULL, tibble_metadata = NULL, replicate_order = NULL, plot_shapes = TRUE, centroid_labels = TRUE, point_size = 5, centroid_label_size = 4, legend_position = "none", PCA_depths_y = NULL, PCA_depths_x = NULL, input_colour_limits = NULL, input_colour_value = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width = 10, height = 10) {
  
  # DEBUG ###
  # table_matrix <- tibble_matrix_absolute_psi_in_sample_replicate_format_with_na
  # condition_order <- temp_condition_names
  # replicate_order <- c("r1", "r2")
  # plot_shapes <- FALSE
  # legend_position <- "none"
  # centroid_labels <- TRUE
  # point_size <- 1
  # centroid_label_size <- 1
  # PCA_depths_y <- c(2, 3, 4)
  # PCA_depths_x <- c(1, 2, 3)
  # input_colour_limits <- tibble_sharp_cluster_mapping$condition_names
  # tibble_metadata <- tibble_sharp_cluster_mapping %>% dplyr::ungroup() %>% dplyr::rename("som_cluster" = "cluster", "condition_names" = "condition_names") %>% dplyr::select(som_cluster, condition_names, Description, Biological_source, Type) %>% dplyr::mutate_at(.vars = c("Description", "Biological_source", "Type"), .funs = function(x) {gsub(x = x, pattern = "[^a-zA-Z0-9 -/]", replacement = "")}) 
  # input_colour_value = rainbow(tibble_sharp_cluster_mapping$cluster %>% max) %>% .[tibble_sharp_cluster_mapping$cluster]
  # save_dir <- R_processing_results_dir
  # save_name <- paste("atlas_totalrna_psisigma_consensus_som_clusterone_with_na_plotly", sep = "")
  # graph_title = "Total RNA/PSI-Sigma w/ replicates, coloured by consensus SOM clusters"
  # width <- 40
  # height <- 40
  ###########
  
  table_matrix <- table_matrix %>% t
  # rownames(table_matrix) <- colnames(table_matrix)
  
  tibble_umap_result <- umap::umap(table_matrix) %>% .$layout %>% 
    tibble::as_tibble(rownames = "condition|replicate", .name_repair = "unique") %>%
    setNames(nm = c("condition|replicate", "V1", "V2")) %>%
    tibble::add_column(
      "condition_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\1"),
      "replicate_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\2")
    ) %>% 
    dplyr::inner_join(., tibble_metadata, by = intersect(colnames(tibble_metadata), colnames(.))[intersect(colnames(tibble_metadata), colnames(.)) %in% c("condition_names", "replicate_names")])
  
  # centroid labels: make a separate file containing the centroid locations
  tibble_centroid_locations <- tibble_umap_result %>% 
    dplyr::group_by(condition_names) %>% 
    dplyr::summarise("centroid_x" = mean(V1),
                     "centroid_y" = mean(V2)) %>%
    add_column("cluster_number" = paste(1:nrow(.)), .after = "condition_names")
  
  # append cluster info
  tibble_umap_result_plot <- dplyr::left_join(tibble_umap_result, tibble_centroid_locations, by = "condition_names")
  
  if (is.null(input_colour_limits) == TRUE) {
    
    input_colour_limits <- c(tibble_umap_result_plot$condition_names %>% unique)
    
  }
  
  if (is.null(input_colour_value) == TRUE) {
    
    input_colour_value <- rainbow(n = (condition_order %>% length))
    
  }
  
  ggplot_plot <- ggplot() +
    (if (plot_shapes == TRUE) {geom_point(aes(x = tibble_umap_result_plot$V1, y = tibble_umap_result_plot$V2, shape = tibble_umap_result_plot$replicate_names, color = tibble_umap_result_plot$condition_names, fill = tibble_umap_result_plot$condition_names), size = point_size) } else {geom_blank(aes(x = tibble_umap_result_plot$V1, y = tibble_umap_result_plot$V2))}) +
    scale_shape_manual(name = "Replicate", values = 1:length(replicate_order)) +
    (if (centroid_labels == TRUE) {geom_text(data = tibble_umap_result_plot, mapping = do.call(what = aes, args = list("x" = tibble_umap_result_plot$centroid_x, "y" = tibble_umap_result_plot$centroid_y, "label" = tibble_umap_result_plot$cluster_number, "color" = tibble_umap_result_plot$condition_names) %>% purrr::splice(tibble_umap_result %>% dplyr::select(-`condition|replicate`, -V1, -V2, -contains("condition_names"), -contains("replicate_names")) %>% purrr::array_tree(margin = 2))), size = centroid_label_size)} else {geom_blank(aes(x = tibble_umap_result_plot$centroid_x, y = tibble_umap_result_plot$centroid_y))}) +
    scale_fill_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    scale_colour_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    # scale_color_brewer(name = "Timepoint", palette = "Spectral", breaks = condition_order, limits = condition_order) +
    ggtitle(paste("UMAP projection\n", graph_title, sep = "")) +
    guides(som_cluster = guide_legend(order = 1)) +
    # guides(size = FALSE) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme_bw() +
    theme(text = element_text(family = "Helvetica"), legend.position = legend_position)
  
  plotly::ggplotly(
    p = ggplot_plot,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = "Timepoint",
    originalData = TRUE,
  ) %>% 
    htmlwidgets::saveWidget(paste(save_dir, "/", save_name, ".html", sep = ""))
  
  ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width, height = height, units = "cm")
  # ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".svg", sep = ""), device = "svg", dpi = 600, width = width, height = height, units = "cm")
  
  write.table(x = tibble_centroid_locations, file = paste(save_dir, "/", save_name, "_cluster_legend.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

## old version
plot_UMAP_for_timepoint_and_replicate_plotly <- function(transposed_matrixtable, timepoint_order = NULL, tibble_colour_metadata = NULL, replicate_order = NULL, plot_shapes = TRUE, centroid_labels = TRUE, point_size = 5, centroid_label_size = 4, legend_position = "none", PCA_depths_y = NULL, PCA_depths_x = NULL, input_colour_limits = NULL, input_colour_value = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width = 10, height = 10, ...) {
  
  # DEBUG ###
  # transposed_matrixtable <- tibble_matrix_absolute_psi_in_sample_replicate_format_with_na %>% t %>% .[gsub(x = rownames(.), pattern = "^(.*)\\|.*$", replacement = "\\1") %in% tibble_sharp_cluster_mapping$condition_names, ]
  # timepoint_order <- temp_condition_names
  # replicate_order <- c("r1", "r2")
  # plot_shapes <- FALSE
  # legend_position <- "none"
  # centroid_labels <- TRUE
  # point_size <- 1
  # centroid_label_size <- 1
  # PCA_depths_y <- c(2, 3, 4)
  # PCA_depths_x <- c(1, 2, 3)
  # input_colour_limits <- tibble_sharp_cluster_mapping$condition_names
  # tibble_colour_metadata <- tibble_sharp_cluster_mapping %>% dplyr::ungroup() %>% dplyr::rename("som_cluster" = "cluster", "input_colour_limits" = "condition_names") %>% dplyr::select(id, som_cluster, input_colour_limits, Description, Biological_source, Type) %>% dplyr::mutate_at(.vars = c("Description", "Biological_source", "Type"), .funs = function(x) {gsub(x = x, pattern = "[^a-zA-Z0-9 -/]", replacement = "")})
  # input_colour_value = rainbow(tibble_sharp_cluster_mapping$cluster %>% max) %>% .[tibble_sharp_cluster_mapping$cluster]
  # save_dir <- R_processing_results_dir
  # save_name <- paste("atlas_polya_psisigma_consensus_som_clusterone_with_na_plotly", sep = "")
  # graph_title = "PolyA/PSI-Sigma w/ replicates, coloured by consensus SOM clusters"
  # width <- 40
  # height <- 40
  ###########
  
  tibble_umap_result <- umap::umap(transposed_matrixtable) %>% .$layout %>% 
    as_tibble(rownames = "condition|replicate", .name_repair = "unique") %>%
    setNames(nm = c("condition|replicate", "V1", "V2")) %>%
    add_column("condition_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\1"),
               "replicate_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\2")
    )
  
  # centroid labels: make a separate file containing the centroid locations
  tibble_centroid_locations <- tibble_umap_result %>% 
    dplyr::group_by(condition_names) %>% 
    dplyr::summarise("centroid_x" = mean(V1),
                     "centroid_y" = mean(V2)) %>%
    add_column("cluster_number" = paste(1:nrow(.)), .after = "condition_names")
  
  tibble_colour_metadata <- tibble_colour_metadata[match(tibble_centroid_locations$condition_names, tibble_colour_metadata$input_colour_limits), ]
  
  # append cluster info
  tibble_umap_result_plot <- dplyr::left_join(tibble_umap_result, tibble_centroid_locations, by = "condition_names")
  
  if (is.null(input_colour_limits) == TRUE) {
    
    input_colour_limits <- c(tibble_umap_result_plot$condition_names %>% unique)
    
  }
  
  if (is.null(input_colour_value) == TRUE) {
    
    input_colour_value <- rainbow(n = (timepoint_order %>% length))
    
  }
  
  # plot UMAP for multiple depths all in one go.
  ggplot_plot <- ggplot() +
    (if (plot_shapes == TRUE) {geom_point(aes(x = tibble_umap_result_plot$V1, y = tibble_umap_result_plot$V2, shape = tibble_umap_result_plot$replicate_names, color = tibble_umap_result_plot$condition_names, fill = tibble_umap_result_plot$condition_names), size = point_size) } else {geom_blank(aes(x = tibble_umap_result_plot$V1, y = tibble_umap_result_plot$V2))}) +
    scale_shape_manual(name = "Replicate", values = 15+(1:length(replicate_order))) +
    (if (centroid_labels == TRUE) {geom_text(data = tibble_centroid_locations, aes(x = centroid_x, y = centroid_y, label = cluster_number, color = condition_names, som_cluster = tibble_colour_metadata$som_cluster, Description = tibble_colour_metadata$Description, Biological_source = tibble_colour_metadata$Biological_source, Type = tibble_colour_metadata$Type), size = centroid_label_size)} else {geom_blank(aes(x = tibble_centroid_locations$centroid_x, y = tibble_centroid_locations$centroid_y))}) +
    scale_fill_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    scale_colour_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    # scale_color_brewer(name = "Timepoint", palette = "Spectral", breaks = timepoint_order, limits = timepoint_order) +
    ggtitle(paste("UMAP projection\n", graph_title, sep = "")) +
    guides(som_cluster = guide_legend(order = 1)) +
    # guides(size = FALSE) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme_bw() +
    theme(text = element_text(family = "Helvetica"), legend.position = legend_position)
  
  plotly::ggplotly(
    p = ggplot_plot,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = "Timepoint",
    originalData = TRUE,
  ) %>% 
    htmlwidgets::saveWidget(paste(save_dir, "/", save_name, ".html", sep = ""))
  
  ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width, height = height, units = "cm")
  # ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".svg", sep = ""), device = "svg", dpi = 600, width = width, height = height, units = "cm")
  
  write.table(x = tibble_centroid_locations, file = paste(save_dir, "/", save_name, "_cluster_legend.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

# FUNCTION TO PLOT tSNE FOR SAMPLE AND REPLICATE
# tibble_metadata: this should have n columns equal to the number of annotation to be specified, and m rows equal to number of samples (including replicates). basically a collection of vectors of length equal to number of samples.
# protected colname: "condition_names", "replicate_names"
# each row is a SAMPLE, each column is a feature
# provide table_matrix as a df with rownames named after each sample
plot_tSNE_condition_and_replicate_plotly <- function(table_matrix, condition_order = NULL, replicate_order = NULL, tibble_metadata = NULL, plot_shapes = TRUE, centroid_labels = TRUE, point_size = 5, centroid_label_size = 4, legend_position = "none", PCA_depths_y = NULL, PCA_depths_x = NULL, input_colour_limits = NULL, input_colour_value = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width = 10, height = 10, .seed = 8, .print_plot = FALSE, .perplexity = NULL, .theta = NULL) {
  
  # DEBUG ###
  # table_matrix = tibble_matrix_absolute_psi_in_sample_replicate_format %>% na.omit %>% t %>% .[, gsub(x = rownames(.), pattern = "^(.*)\\|.*$", replacement = "\\1") %in% tibble_sharp_cluster_mapping$condition_names]
  # condition_order = temp_condition_names
  # replicate_order = c("r1", "r2")
  # tibble_metadata <- tibble_sharp_cluster_mapping %>% dplyr::ungroup() %>% dplyr::rename("som_cluster" = "cluster", "condition_names" = "condition_names") %>% dplyr::select(som_cluster, condition_names, Description, Biological_source, Type) %>% dplyr::mutate_at(.vars = c("Description", "Biological_source", "Type"), .funs = function(x) {gsub(x = x, pattern = "[^a-zA-Z0-9 -/]", replacement = "", useBytes = TRUE)})
  # plot_shapes = FALSE
  # legend_position = "none"
  # centroid_labels = TRUE
  # point_size = 1
  # centroid_label_size = 1
  # PCA_depths_y = c(2, 3, 4)
  # PCA_depths_x = c(1, 2, 3)
  # input_colour_limits = tibble_sharp_cluster_mapping$condition_names
  # input_colour_value = rainbow(tibble_sharp_cluster_mapping$cluster %>% max) %>% .[tibble_sharp_cluster_mapping$cluster]
  # save_dir = R_processing_results_dir
  # save_name = paste(vector_experiment_metadata_initial %>% paste(collapse = "_"), "_tsne_consensus_som_clusterone_replicates_no_na_", sampletype, sep = "")
  # graph_title = "Total RNA/Leafcutter w/ replicates, coloured by consensus SOM clusters"
  # width = 40
  # height = 40
  ###########
  
  set.seed(.seed)
  
  df_tsne_result <- callr_insulator(.f = quote( 
    Rtsne::Rtsne(X = table_matrix %>% unique, check_duplicates = FALSE, verbose = TRUE, num_threads = 0, perplexity = .perplexity, theta = .theta) 
  ) ) %>% 
    .$Y %>% 
    as.data.frame() %>%
    dplyr::mutate("condition|replicate" = rownames(table_matrix %>% unique), .before = 1) %>%
    dplyr::mutate(
      "condition_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\1"),
      "replicate_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\2")
    ) %>% 
    dplyr::inner_join(., tibble_metadata, by = intersect(colnames(tibble_metadata), colnames(.))[intersect(colnames(tibble_metadata), colnames(.)) %in% c("condition_names", "replicate_names")])
  
  # centroid labels: make a separate file containing the centroid locations
  tibble_centroid_locations <- df_tsne_result %>% 
    dplyr::group_by(condition_names) %>% 
    dplyr::summarise("centroid_x" = mean(V1),
                     "centroid_y" = mean(V2)) %>%
    dplyr::mutate("cluster_number" = paste(1:nrow(.)), .after = "condition_names")
  
  # append cluster info
  df_tsne_result_plot <- dplyr::left_join(df_tsne_result, tibble_centroid_locations, by = "condition_names")
  
  if (is.null(input_colour_limits) == TRUE) {
    input_colour_limits <- c(df_tsne_result_plot$condition_names %>% unique)
  }
  
  if (is.null(input_colour_value) == TRUE) {
    input_colour_value <- rainbow(n = (timepoint_order %>% length))
  }
  
  ggplot_plot <- ggplot() +
    (if (plot_shapes == TRUE) {geom_point(aes(x = df_tsne_result_plot$V1, y = df_tsne_result_plot$V2, shape = df_tsne_result_plot$replicate_names, color = df_tsne_result_plot$condition_names, fill = df_tsne_result_plot$condition_names), size = point_size) } else {geom_blank(aes(x = df_tsne_result_plot$V1, y = df_tsne_result_plot$V2))}) +
    scale_shape_manual(name = "Replicate", values = 1:length(replicate_order)) +
    (if (centroid_labels == TRUE) {geom_text(data = df_tsne_result_plot, mapping = do.call(what = aes, args = list("x" = df_tsne_result_plot$centroid_x, "y" = df_tsne_result_plot$centroid_y, "label" = df_tsne_result_plot$cluster_number, "color" = df_tsne_result_plot$condition_names) %>% purrr::splice(df_tsne_result_plot %>% dplyr::select(-`condition|replicate`, -V1, -V2, -contains("condition_names"), -contains("replicate_names")) %>% purrr::array_tree(margin = 2))), size = centroid_label_size)} else {geom_blank(aes(x = df_tsne_result_plot$centroid_x, y = df_tsne_result_plot$centroid_y))}) +
    scale_fill_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    scale_colour_manual(name = "Timepoint", breaks = input_colour_limits, limits = input_colour_limits, values = input_colour_value) +
    # scale_color_brewer(name = "Timepoint", palette = "Spectral", breaks = timepoint_order, limits = timepoint_order) +
    ggtitle(paste("UMAP projection\n", graph_title, sep = "")) +
    guides(som_cluster = guide_legend(order = 1)) +
    # guides(size = FALSE) +
    xlab("tSNE_1") +
    ylab("tSNE_2") +
    theme_bw() +
    theme(text = element_text(family = "Helvetica"), legend.position = legend_position)
  
  if (.print_plot == TRUE) {
    print(ggplot_plot)
  }
  
  plotly::ggplotly(
    p = ggplot_plot,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = "Timepoint",
    originalData = TRUE,
  ) %>% 
    htmlwidgets::saveWidget(paste(save_dir, "/", save_name, ".html", sep = ""))
  
  ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width, height = height, units = "cm")
  # ggsave(plot = ggplot_plot, filename = paste(save_dir, "/", save_name, ".svg", sep = ""), device = "svg", dpi = 600, width = width, height = height, units = "cm")
  
  write.table(x = tibble_centroid_locations, file = paste(save_dir, "/", save_name, "_cluster_legend.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

plot_PCA_for_sample_and_replicate <- function(matrixtable, timepoint_order = NULL, replicate_order = NULL, plot_shapes = TRUE, text_labels = FALSE, point_or_label_size = 2, legend_position = "right", PCA_depths_y = NULL, PCA_depths_x = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width = 10, height = 10) {
  
  # DEBUG ###
  # matrixtable <- wide_tibble_matrix_processed_sorted_tibbles %>% as.data.frame
  # save_dir <- R_processing_results_dir
  # save_name <- "total_RNA"  # 
  # graph_title <- "Total RNA"  # 
  # timepoint_order <- long_tibble_processed_sorted_tibbles$sample_name %>% unique
  # replicate_order <- c("r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8")
  # width = 10
  # height = 10
  # PCA_depths_y = c(2, 3, 4)
  # PCA_depths_x = c(1, 2, 3)
  # text_labels <- FALSE
  # point_or_label_size <- 4
  ###########
  
  PCA_result <- prcomp(matrixtable)
  
  ## measure PC variance #
  PCA_stdev <- tibble("PC" = 1:(PCA_result[["sdev"]] %>% length), "stdev" = PCA_result[["sdev"]])
  
  PCA_variance <- tibble(PC = PCA_stdev$PC, variance = PCA_stdev$stdev ^ 2)
  PCA_variance <- add_column(PCA_variance, variance_explained = PCA_variance$variance/sum(PCA_variance$variance) * 100)
  
  ggplot(PCA_variance) + 
    geom_col(aes(x = PC, y = variance_explained, fill = PC)) +
    scale_fill_gradientn(colours = heat.colors(n = (PCA_result[["sdev"]] %>% length))) +
    ggtitle(paste("PCA variance distribution\n", graph_title, sep = "")) +
    guides(size = FALSE) + 
    xlab("PC") +
    ylab("Variance explained (%)") +
    theme_bw() +
    theme(text = element_text(family = "Helvetica")) +
    ggsave(filename = paste(save_dir, "PCA_barplot_stdevs_", save_name, ".pdf", sep = ""), device = "pdf", dpi = 600, width = 45, height = 25, units = "cm") +
    ggsave(filename = paste(save_dir, "PCA_barplot_stdevs_", save_name, ".svg", sep = ""), device = "svg", dpi = 600, width = 45, height = 25, units = "cm")
  
  ## measure PC loadings #
  ## column names of the matrix need to be split by a "|", like this: timepoint|replicatenumber
  PCA_loadings <- PCA_result[["rotation"]] %>% as_tibble(rownames = "sample")
  PCA_loadings <- add_column(PCA_loadings, "timepoint" = gsub(x = PCA_loadings$sample, pattern = "(.*)\\|(.*)", replacement = "\\1") %>% factor(x = ., levels = timepoint_order), .after = "sample")
  PCA_loadings <- add_column(PCA_loadings, "replicatenumber" = gsub(x = PCA_loadings$sample, pattern = "(.*)\\|(.*)", replacement = "\\2") %>% factor(x = ., levels = replicate_order), .after = "sample")
  PCA_loadings <- add_column(PCA_loadings, 
                             "condition_names" = gsub(x = PCA_loadings$sample, pattern = "(.*)\\|(.*)", replacement = "\\1"),
                             "replicate_names" = gsub(x = PCA_loadings$sample, pattern = "(.*)\\|(.*)", replacement = "\\2"), 
                             .after = "sample")
  # append cluster number according to the order of timepoints specified by the user
  PCA_loadings <- dplyr::left_join(PCA_loadings, tibble("condition_names" = levels(PCA_loadings$timepoint), "cluster_number" = 1:length(levels(PCA_loadings$timepoint))), by = "condition_names") %>% 
    dplyr::relocate(cluster_number, .after = "sample")
  
  # plot PCA for multiple depths all in one go.
  purrr::map2(.x = PCA_depths_x, .y = PCA_depths_y, .f = function(.x, .y) {
    
    pc_x <- .x
    pc_y <- .y
    
    ggplot(PCA_loadings) + 
      (if (text_labels == FALSE) {geom_point(aes(x = !!(paste("PC", pc_x, sep = "") %>% as.name), y = !!(paste("PC", pc_y, sep = "") %>% as.name), shape = replicatenumber, color = timepoint), size = point_or_label_size)} else if (text_labels == TRUE) {geom_text(aes(x = !!(paste("PC", pc_x, sep = "") %>% as.name), y = !!(paste("PC", pc_y, sep = "") %>% as.name), label = PCA_loadings$cluster_number, color = timepoint), size = point_or_label_size, position = position_dodge(width = 4))}) +
      scale_color_brewer(name = "Timepoint", palette = "Spectral", breaks = timepoint_order, limits = timepoint_order) +
      scale_shape_manual(name = "Replicate", values = 1:length(replicate_order)) +
      ggtitle(paste("PCA loadings\n", graph_title, sep = "")) +
      guides(size = FALSE) + 
      xlab(paste("PC", pc_x, " (", PCA_variance[pc_x, "variance_explained"] %>% signif(3), "%)", sep = "")) +
      ylab(paste("PC", pc_y, " (", PCA_variance[pc_y, "variance_explained"] %>% signif(3), "%)", sep = "")) +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"), legend.position = legend_position) +
      ggsave(filename = paste(save_dir, "PCA_loadings_", save_name, "_PC", pc_y, "_vs_PC_", pc_x, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width, height = height, units = "cm") +
      ggsave(filename = paste(save_dir, "PCA_loadings_", save_name, "_PC", pc_y, "_vs_PC_", pc_x, ".svg", sep = ""), device = "svg", dpi = 600, width = width, height = height, units = "cm")
    
  } )
  
}

# new, revamped version
## column names of the table_matrix need to be split by a "|", like this: condition|replicate
## each column is a sample and each row is a feature
## tibble_metadata must contain no. rows the same number as the number of columns in the input data
### protected column names: ".sample" in the format condition|replicate, each row specifies additional annotations for each data point
ggplot_plotly_pca <- function(table_matrix, .scale = TRUE, condition_order = NULL, replicate_order = NULL, tibble_metadata = NULL, plot_shapes = TRUE, centroid_labels = TRUE, point_size = 5, centroid_label_size = 4, legend_position = "none", pca_depths_x = NULL, pca_depths_y = NULL, ggplot_colour_limits = NULL, ggplot_colour_values = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width_embeddings_plot = 10, height_embeddings_plot = 10, width_variance_plot = 45, height_variance_plot = 25, .plot_inline = TRUE, .save_plots = TRUE) {
  
  # DEBUG ###
  # table_matrix <- matrix_donor_snp_presence[1:1000, ] %>% (function (x) {colnames(x) <- temp_tibble_metadata$`.sample`; return(x)})
  # .scale <- TRUE
  # condition_order <- tibble_donor_pair_record_processed$donor_batch %>% unique
  # replicate_order <- c(1, 2, 3)
  # tibble_metadata <- temp_tibble_metadata
  # plot_shapes <- TRUE
  # centroid_labels <- TRUE
  # point_size <- 5
  # centroid_label_size <- 4
  # legend_position <- "none"
  # pca_depths_x <- c(1, 2, 3)
  # pca_depths_y <- c(2, 3, 4)
  # ggplot_colour_limits <- NULL
  # ggplot_colour_values <- NULL
  # save_dir <- r_results_dir
  # save_name <- c(vector_experiment_tag, "donor_snp_presence_pca_plot") %>% paste(collapse = "_")
  # graph_title <- NULL
  # width_embeddings_plot <- 10
  # height_embeddings_plot <- 10
  # width_variance_plot <- 45
  # height_variance_plot <- 25
  # .plot_inline <- TRUE
  # .save_plots <- TRUE
  # ###########
  
  # result: stdevs + a frame of rotations (embeddings) with columns PCs and rows as each sample
  # pca_result <- prcomp(table_matrix, scale. = .scale, rank. = max(c(pca_depths_y, pca_depths_x)))
  # summary.prcomp_result <- summary(pca_result)
  
  # irlba because it scales well
  irlba_prcomp_result <- irlba::prcomp_irlba(x = table_matrix, n = max(c(pca_depths_y, pca_depths_x)), retx = FALSE, scale. = .scale)
  vector_variances_explained <- ((irlba_prcomp_result$sdev)^2)/sum((irlba_prcomp_result$sdev)^2)
  
  # get PC variance
  # first get stdev
  # tibble_pca_stdev <- tibble::tibble("PC" = 1:(pca_result[["sdev"]] %>% length), "stdev" = pca_result[["sdev"]])
  # now get variance
  # tibble_pca_variance <- tibble::tibble(PC = tibble_pca_stdev$PC, variance = tibble_pca_stdev$stdev ^ 2)
  # tibble_pca_variance <- tibble::add_column(tibble_pca_variance, variance_explained = tibble_pca_variance$variance/sum(tibble_pca_variance$variance) * 100)
  
  # plot PC variance
  # ggplot_variance <- ggplot2::ggplot(tibble_pca_variance) + 
  # ggplot2::geom_col(aes(x = PC, y = variance_explained, fill = PC)) +
  # ggplot2::scale_fill_gradientn(colours = heat.colors(n = (pca_result[["sdev"]] %>% length))) +
  # ggplot2::ggtitle(paste("pca variance distribution\n", graph_title, sep = "")) +
  # ggplot2::guides(size = FALSE) + 
  # ggplot2::xlab("PC") +
  # ggplot2::ylab("Variance explained (%)") +
  # ggplot2::theme_bw() +
  # ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))
  
  # if (.save_plots == TRUE) {
  #   ggplot2::ggsave(plot = ggplot_variance, filename = paste(save_dir, "pca_barplot_stdevs_", save_name, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width_variance_plot, height = height_variance_plot, units = "cm")
  # }
  
  # if (.plot_inline == TRUE) {
  #   print(ggplot_variance)
  # }
  
  # get PC loadings
  ## column names of the matrix need to be split by a "|", like this: condition|replicate
  # tibble_pca_loadings <- pca_result[["rotation"]] %>% tibble::as_tibble(rownames = ".sample")
  tibble_pca_loadings <- irlba_prcomp_result[["rotation"]] %>% tibble::as_tibble() %>% dplyr::mutate(".sample" = tibble_metadata$.sample)
  tibble_pca_loadings <- dplyr::mutate(
    tibble_pca_loadings, 
    "condition" = gsub(x = tibble_pca_loadings$`.sample`, pattern = "(.*)\\|(.*)", replacement = "\\1") %>% factor(x = ., levels = condition_order),
    "replicate" = gsub(x = tibble_pca_loadings$`.sample`, pattern = "(.*)\\|(.*)", replacement = "\\2") %>% factor(x = ., levels = replicate_order),
    .after = ".sample")
  # add tibble_metadata
  tibble_pca_loadings <- dplyr::left_join(
    tibble_pca_loadings, 
    tibble_metadata, 
    by = ".sample"
  )
  
  # automatically make color scale if not specified by the user
  if (is.null(ggplot_colour_limits) == TRUE) {
    ggplot_colour_limits <- levels(tibble_pca_loadings$condition)
  }
  
  if (is.null(ggplot_colour_values) == TRUE) {
    ggplot_colour_values <- rainbow(n = length(levels(tibble_pca_loadings$condition)))
  }
  
  # plot pca for multiple depths all in one go.
  purrr::map2(
    .x = pca_depths_x, 
    .y = pca_depths_y, 
    .f = function(a1, a2) {
      
      # DEBUG ###
      # a1 <- 1
      # a2 <- 2
      ###########
      
      pc_x <- a1
      pc_y <- a2
      
      # centroid labels: make a separate table containing the centroid locations
      tibble_centroid_locations <- tibble_pca_loadings %>% 
        dplyr::group_by(condition) %>% 
        dplyr::summarise("centroid_x" = mean(!!as.name(paste("PC", pc_x, sep = ""))),
                         "centroid_y" = mean(!!as.name(paste("PC", pc_y, sep = "")))) %>%
        dplyr::mutate("centroid_number" = paste(1:nrow(.)), .after = "condition")
      
      # join centroid table onto the main table
      tibble_pca_loadings_with_centroids <- dplyr::left_join(
        tibble_pca_loadings,
        tibble_centroid_locations,
        by = "condition"
      )
      
      ggplot_embedding <- ggplot2::ggplot() + 
        # shapes
        (if (plot_shapes == TRUE) {ggplot2::geom_point(mapping = do.call(what = ggplot2::aes, args = list("x" = tibble_pca_loadings_with_centroids$centroid_x, "y" = tibble_pca_loadings_with_centroids$centroid_y, "shape" = tibble_pca_loadings_with_centroids$replicate, "color" = tibble_pca_loadings_with_centroids$condition, "fill" = tibble_pca_loadings_with_centroids$condition) %>% purrr::splice(tibble_pca_loadings_with_centroids %>% dplyr::select(-`.sample`, -contains(match = "PC", ignore.case = FALSE), -condition, -replicate) %>% purrr::array_tree(margin = 2))), size = point_size) } else {ggplot2::geom_blank(aes(x = tibble_pca_loadings_with_centroids$V1, y = tibble_pca_loadings_with_centroids$V2))}) +
        scale_shape_manual(name = "Replicate", values = 1:length(replicate_order)) +
        # centroid labels
        (if (centroid_labels == TRUE) {ggplot2::geom_text(data = tibble_pca_loadings_with_centroids, mapping = do.call(what = ggplot2::aes, args = list("x" = tibble_pca_loadings_with_centroids$centroid_x, "y" = tibble_pca_loadings_with_centroids$centroid_y, "label" = tibble_pca_loadings_with_centroids$centroid_number, "color" = tibble_pca_loadings_with_centroids$condition) %>% purrr::splice(tibble_pca_loadings_with_centroids %>% dplyr::select(-`.sample`, -contains(match = "PC", ignore.case = FALSE), -condition, -replicate) %>% purrr::array_tree(margin = 2))), size = centroid_label_size)} else {geom_blank(aes(x = tibble_pca_loadings_with_centroids$centroid_x, y = tibble_pca_loadings_with_centroids$centroid_y))}) +
        
        scale_fill_manual(name = "Condition", breaks = ggplot_colour_limits, limits = ggplot_colour_limits, values = ggplot_colour_values) +
        scale_colour_manual(name = "Condition", breaks = ggplot_colour_limits, limits = ggplot_colour_limits, values = ggplot_colour_values) +
        
        ggtitle(paste("PCA loadings\n", graph_title, sep = "")) +
        
        # xlab(paste("PC", pc_x, " (", 100 * summary.prcomp_result$importance %>% .[2, pc_x] %>% signif(3), "%)", sep = "")) +
        # ylab(paste("PC", pc_y, " (", 100 * summary.prcomp_result$importance %>% .[2, pc_y] %>% signif(3), "%)", sep = "")) +
        
        xlab(paste("PC", pc_x, " (", 100 * vector_variances_explained %>% .[pc_x] %>% signif(3), "%)", sep = "")) +
        ylab(paste("PC", pc_y, " (", 100 * vector_variances_explained %>% .[pc_y] %>% signif(3), "%)", sep = "")) +
        
        theme_bw() +
        theme(text = element_text(family = "Helvetica"), legend.position = legend_position)
      
      if (.plot_inline == TRUE) {
        print(ggplot_embedding)
      }
      
      if (.save_plots == TRUE) {
        
        ggsave(plot = ggplot_embedding, filename = paste(save_dir, "/", save_name, "_pc", pc_y, "_vs_pc_", pc_x, ".pdf", sep = ""), device = "pdf", dpi = 600, width = width_embeddings_plot, height = height_embeddings_plot, units = "cm")
        
        plotly::ggplotly(
          p = ggplot_embedding,
          width = NULL,
          height = NULL,
          tooltip = "all",
          dynamicTicks = FALSE,
          layerData = "condition",
          originalData = TRUE,
        ) %>% 
          htmlwidgets::saveWidget(paste(save_dir, "/", save_name, "_pc", pc_y, "_vs_pc_", pc_x, ".html", sep = ""))
        
        write.table(x = tibble_pca_loadings_with_centroids, file = paste(save_dir, "/", save_name, "_legend.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        
      }
      
    } )
  
}

# sequential splitter by segment size
## split list/vector x into whatever many segments of size n
## keep_remainders: whether to keep the last partially filled element due to incomplete division
seqsplitter <- function(x, n, keep_remainders = TRUE) {
  
  if (keep_remainders == TRUE) {
    list_result <- purrr::map(
      .x = (1:((length(x)%/%n) + 1)), 
      .f = ~x[(n*(.x - 1) + 1):min(c(length(x), n*.x))])
  } else {
    list_result <- purrr::map(
      .x = 1:(length(x)%/%n), 
      .f = ~x[(n*(.x - 1) + 1):min(c(length(x), n*.x))])
  }
  
  return(list_result)
  
}

## takes a tibble of x and y coords along with a grouping label. colnames: x, y, group.
## number_of_output_datapoints: number of points to specify in the output curve

## The function will automatically chop up the max/min x-values and loop through each of these output datapoints. 
## for each datapoint, we loop through each group
## for each group, compute the linearly interpolated y-value based on the closest value to the left/right of the datapoint.
## Finally, average the interpolated value for each group
sum_multiple_curves_from_datapoints <- function(input_tibble = NULL, number_of_output_datapoints = NULL, ncores_L1 = 2, ncores_L2 = 1) {
  
  # DEBUG ###
  # input_tibble <- tibble(
  # "y" = purrr::map(.x = multiclass.roc_result$rocs, .f = ~.x$sensitivities) %>% unlist, 
  # "x" = purrr::map(.x = multiclass.roc_result$rocs, .f = ~(1 - .x$specificities)) %>% unlist, 
  # "group" = purrr::imap(.x = multiclass.roc_result$rocs , .f = ~rep(.y, times = length(.x$specificities))) %>% unlist %>% paste("g", ., sep = "")) 
  # number_of_output_datapoints <- 100
  # ncores_L1 <- 4
  # ncores_L2 <- 1
  ###########
  
  # split tibble by group
  input_tibble_split_by_group <- input_tibble %>% dplyr::group_split(group) %>% purrr::set_names(nm = purrr::map(.x = ., .f = ~.x$group %>% unique) %>% unlist)
  
  # split interval into desired output datapoints
  vector_output_x_coords <- seq(from = min(input_tibble$x), to = max(input_tibble$y), length.out = number_of_output_datapoints)
  
  plan(list(tweak(multicore, workers = ncores_L1),
            tweak(multicore, workers = ncores_L2)))
  
  # loop through each output datapoint and linearly interpolate
  list_averaged_y_values_at_each_output_datapoint <- furrr::future_map(
    .x = vector_output_x_coords,
    .f = function(a1) {
      
      # DEBUG ###
      # a1 <- vector_output_x_coords[[1]]
      ###########
      
      list_interpolated_y_values_per_group <- furrr::future_map(
        .x = input_tibble_split_by_group,
        .f = function(b1) {
          
          approx(x = b1$x, y = b1$y, xout = a1, n = 1, ties = mean) %>%
            .$y %>% 
            return
          
        } )
      
      # return the average of interpolated y values
      return(
        mean(list_interpolated_y_values_per_group %>% unlist)
      )
      
    }, .progress = TRUE )
  
  # return tibble of averaged values
  return(
    tibble(
      "x" = vector_output_x_coords,
      "y" = list_averaged_y_values_at_each_output_datapoint %>% unlist
    )
  )
  
}
