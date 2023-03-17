# This document contains all the amazing scripts necessary to do regular bioinformatics manipulations.

# simple round robin callr
## true round robin
## one line per process
## global environment automatically exported via temp file
round_robin_pmap_callr <- function(.l, .f, .num_workers = 1, .temp_path = NULL, .re_export = TRUE, .status_messages_dir, .job_name, ...) {
  
  # DEBUG ###
  # .l <- list("a1" = 1:100)
  # .f <- function(a1) {Sys.sleep(30); return(a1 + 1)}
  # # .f <- function(a1, a2) {return(a1 + a2)}
  # .num_workers <- 8
  # .temp_path <- "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/temp/atlas_polya_psisigma_exons_3FT_unmatched_list_131_temp.RData"
  # .status_messages_dir <- "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/temp/"
  # .job_name <- output_name
  # .re_export <- FALSE
  
  # .l <- list(
  #   "a1" = 1:length(list_alt_exons_sectored),
  #   "a2" = purrr::map2(.x = list_island_intervals, .y = vector_ref_genome_fasta_path, .f = ~rep(.y, times = nrow(.x))) %>% flatten
  #   # setdiff(1:length(list_alt_exons_sectored), list.files(paste(output_dir, "/temp/", sep = ""), pattern = "list_three_frame_translation") %>% gsub(pattern = ".*_(\\d+)_temp.RData", replacement = "\\1", perl = TRUE) %>% type.convert(as.is = TRUE))
  # )
  # .f <- testfun
  # .num_workers <- 8
  # .temp_path <- paste(output_dir, "/temp/", output_name, "_callr_export_3FT_temp.RData", sep = "")
  # .status_messages_dir <- paste(output_dir, "temp/", sep = "")
  # .job_name <- output_name
  # .re_export <- FALSE
  ###########
  
  if (is.null(.job_name)) {
    .job_name <- as.numeric(Sys.time())
  }
  
  print(.temp_path)
  
  print(file.exists(.temp_path))
  print(.re_export)
  
  print((!file.exists(.temp_path)) | .re_export == TRUE)
  
  # save the current env into temp path
  if (!file.exists(.temp_path) | .re_export == TRUE) {
    save.image(file = .temp_path, ...)
  }
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(purrr::map(.x = .l, .f = ~length(.x))))
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  # preallocate worker list
  list_workers <- list()
  
  # initial allocation of tasks to maximum number of workers
  for (..i in 1:.num_workers) {
    
    function_current_instance <- .f
    
    list_current_arguments <- purrr::map(.x = .l, .f = ~.x[[..i]])
    
    formals(fun = function_current_instance) <- list_current_arguments
    
    # essential to spike the exported function with a command to read the environment back in
    body(function_current_instance) <- as.call(purrr::prepend(as.list(body(function_current_instance)), expression(load(file = commandArgs()[2])), before = 2))
    
    assign(
      x = paste("worker_", ..i, sep = ""), 
      value = callr::r_bg(
        cmdargs	= .temp_path,
        func = function_current_instance
      )
    )
    
    list_workers <- purrr::splice(list_workers, get(paste("worker_", ..i, sep = "")))
    names(list_workers)[length(list_workers)] <- paste("worker_", ..i, sep = "")
    
  }
  
  # keep track of iteration progress
  current_map_index <- .num_workers
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  while(all(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0)) {
    
    vector_logical_indices_workers_completed_reported <- unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0
    
    # check on process status and write progress file
    purrr::map2(
      .x = list_workers, 
      .y = names(list_workers), 
      .f = function(a1, a2) {
        
        stdout_lines <- a1$read_output_lines()
        stderr_lines <- a1$read_error_lines()
        
        if (length(stdout_lines) != 0) {write.table(x = stdout_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stdout.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)}
        
        if (length(stderr_lines) != 0) {write.table(x = stderr_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)}
        
      } )
    
    # round robin action until the list is fully mapped
    if (current_map_index < map_length) {
      
      number_of_processes_alive <- length(which(unlist(purrr::map(.x = list_workers, .f = ~.x$is_alive()))))
      
      if (number_of_processes_alive < .num_workers) {
        
        new_map_start <- current_map_index + 1
        new_map_end <- min(c(current_map_index + .num_workers - number_of_processes_alive, map_length))
        
        current_map_index <- new_map_end
        
        for (..i in (new_map_start):min(c(new_map_end, map_length))) {
          
          function_current_instance <- .f
          worker_1$read_output_lines()
          list_current_arguments <- purrr::map(.x = .l, .f = ~.x[[..i]])
          
          formals(fun = function_current_instance) <- list_current_arguments
          
          # essential to spike the exported function with a command to read the environment back in
          body(function_current_instance) <- as.call(prepend(as.list(body(function_current_instance)), expression(load(file = commandArgs()[2])), before = 2))
          
          assign(
            x = paste("worker_", ..i, sep = ""), 
            value = callr::r_bg(
              cmdargs	= .temp_path,
              func = function_current_instance
            )
          )
          
          list_workers <- purrr::splice(list_workers, get(paste("worker_", ..i, sep = "")))
          names(list_workers)[length(list_workers)] <- paste("worker_", ..i, sep = "")
          
        }
        
      }
      
      # an actually useful progress bar although rudimentary
      cat(paste("\rPercent map completion: ", current_map_index, "/", map_length, " (", round(x = 100*current_map_index/map_length, digits = 1), "%)", sep = ""))
      
      Sys.sleep(1)
      
    } else if (current_map_index == map_length) {
      
      if (length(which(vector_logical_indices_workers_completed_reported)) == map_length) {
        break()
      }
      
    }
    
  }
  
  if (any(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)) {
    purrr::map(.x = list_workers, .f = ~.x$signal(9))
    stop(print("Exit status failure received on workers: ", paste(names(list_workers)[which(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)]), collapse = ","))
  }
  
  rm(list_workers)
  rm(list = ls(pattern = "worker_"))
  
  cat("\n")
  
}