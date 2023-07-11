# This document contains all the amazing scripts necessary to do regular bioinformatics manipulations.

# simple round robin callr
## true round robin
## one line per process
## global environment automatically exported via temp file
## three personalities are permitted using the option .result_mode: 1. save - do not return result, just save to disk; 2. ordered - return result but in order (at each tick, check if we have consecutive saved files ready before splicing them into a single list); 3. unordered - return list elements as they arrive

# test_tempfile <- tempfile()
# test_tempdir <- tempdir()
# 
# test <- round_robin_pmap_callr(
#   .l = list(
#     "b1" = a1 %>% purrr::array_tree() %>% head,
#     "b2" = 1:length(a1 %>% purrr::array_tree() %>% head)
#   ),
#   .num_workers = 4,
#   .env_flag = "user",
#   .re_export = TRUE,
#   .temp_path = paste(tempdir, "/tempdata.rdata", sep = ""),
#   .temp_dir = tempdir,
#   .objects = c(),
#   .status_messages_dir = paste(tempdir, sep = ""),
#   .job_name = "test",
#   .result_mode = "unordered",
#   .f = function(b1, b2) {Sys.sleep(0); return(list(LETTERS[b2]))})

round_robin_pmap_callr <- function(.l, .f, .num_workers = 1, .temp_path = NULL, .temp_dir = NULL, .re_export = FALSE, .env_flag = "global", .objects, .status_messages_dir, .job_name, .result_mode = "save", .keep_intermediate_files = FALSE, ...) {
  
  # DEBUG ###
  # .l = list(
  #   "b1" = 1:20
  # )
  # .num_workers = 4
  # .env_flag = "user"
  # .re_export = TRUE
  # .temp_path = paste(tempdir, "/tempdata.rdata", sep = "")
  # .temp_dir = tempdir
  # .objects = c()
  # .status_messages_dir = paste(tempdir, sep = "")
  # .job_name = "test"
  # .result_mode = "unordered"
  # .f = function(b1) {Sys.sleep(120); return(list(LETTERS[b1]))}
  ###########
  
  system("ulimit -n 65536")
  
  if (is.null(.job_name)) {
    .job_name <- as.numeric(Sys.time())
  }
  
  # save the current env into temp path
  if (.re_export == FALSE & file.exists(.temp_path)) {
    message("Temp object data file found. Will load that instead of re-exporting")
  } else {
    message("Exporting object data file")
    
    if (.env_flag == "global") {
      save.image(file = .temp_path)
    } else if (.env_flag == "user") {
      save(list = .objects, file = .temp_path)
    }
    
  }
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(purrr::map(.x = .l, .f = ~length(.x))))
  
  if (.num_workers > map_length) {
    .num_workers <- map_length
  }
  
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
    body(function_current_instance) <- as.call(purrr::prepend(as.list(body(function_current_instance)), str2expression(paste("load(file = \"", .temp_path, "\")", sep = "")), before = 2))
    
    # modify the function depending on what is specified
    if (.result_mode %in% c("ordered", "unordered")) {
      
      function_to_run <- function(f) {
        return(NULL)
      }
      
      body(function_to_run) <- as.call(purrr::prepend(as.list(body(function_to_run)), str2expression(paste("saveRDS(object = f(), file = \"", paste(.temp_dir, "/", .job_name, "_chunk_", ..i, ".rds\")", sep = ""), sep = "")), before = 2))
      
      assign(
        x = paste("chunk_", ..i, sep = ""), 
        value = callr::r_bg(
          cmdargs	= c(.temp_path, ..i),
          args = list("f" = function_current_instance),
          func = function_to_run
        )
      )
      
    } else {
      function_to_run <- function_current_instance
      
      assign(
        x = paste("chunk_", ..i, sep = ""), 
        value = callr::r_bg(
          cmdargs	= c(.temp_path, ..i),
          func = function_to_run
        )
      )
    }
    
    list_workers <- purrr::splice(list_workers, get(paste("chunk_", ..i, sep = "")))
    names(list_workers)[length(list_workers)] <- paste("chunk_", ..i, sep = "")
    
  }
  
  # keep track of iteration progress
  current_map_index <- .num_workers
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  flag_completion <- FALSE
  
  while(all(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0)) {
    
    vector_logical_indices_workers_completed_reported <- unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0
    
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
          write.table(x = stdout_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stdout.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)
          }
        
        if (length(stderr_lines) != 0) {
          write(x = paste(Sys.time()), file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE)
          write.table(x = stderr_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)
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
          
          # print("..i")
          # print(..i)
          
          function_current_instance <- .f
          
          list_current_arguments <- purrr::map(.x = .l, .f = ~.x[[..i]])
          
          formals(fun = function_current_instance) <- list_current_arguments
          
          # essential to spike the exported function with a command to read the environment back in
          body(function_current_instance) <- as.call(purrr::prepend(as.list(body(function_current_instance)), str2expression(paste("load(file = \"", .temp_path, "\")", sep = "")), before = 2))
          
          # modify the function depending on what is specified
          if (.result_mode %in% c("ordered", "unordered")) {
            
            function_to_run <- function(f) {
              return(NULL)
            }
            
            body(function_to_run) <- as.call(purrr::prepend(as.list(body(function_to_run)), str2expression(paste("saveRDS(object = f(), file = \"", paste(.temp_dir, "/", .job_name, "_chunk_", ..i, ".rds\")", sep = ""), sep = "")), before = 2))
            
            assign(
              x = paste("chunk_", ..i, sep = ""), 
              value = callr::r_bg(
                cmdargs	= c(.temp_path, ..i),
                args = list("f" = function_current_instance),
                func = function_to_run
              )
            )
            
          } else {
            function_to_run <- function_current_instance
            
            assign(
              x = paste("chunk_", ..i, sep = ""), 
              value = callr::r_bg(
                cmdargs	= c(.temp_path, ..i),
                func = function_to_run
              )
            )
          }
          
          list_workers <- purrr::splice(list_workers, get(paste("chunk_", ..i, sep = "")))
          names(list_workers)[length(list_workers)] <- paste("chunk_", ..i, sep = "")
          
        }
        
      }
      
      # SPLICER
      
      if (.result_mode %in% c("ordered", "unordered")) {
      
        # initialise list_result if it hasnt already been created
        if (ls(pattern = "^list_result$") %>% length == 0) {
          list_result <- list()
          
          vector_current_chunks_spliced <- 0
        }
        
        if (length(list_result) < map_length) {
          
          vector_completed_chunks <- names(vector_logical_indices_workers_completed_reported) %>% gsub(pattern = "chunk_", replacement = "") %>% type.convert(as.is = TRUE)
          
          vector_chunks_to_be_spliced <- setdiff(vector_completed_chunks, vector_current_chunks_spliced) %>% sort
          
          # consider order
          if (length(vector_chunks_to_be_spliced) > 0) {
            
            # if order is important, simply drop non-consecutive terms from the vector_chunks_to_be_spliced. they will reappear later when the terms in the middle are available.
            if (.result_mode == "ordered") {
              
              global_vector_current_chunks_spliced <<- vector_current_chunks_spliced
              print("vector_current_chunks_spliced")
              print(vector_current_chunks_spliced)
              
              global_vector_completed_chunks <<- vector_completed_chunks
              print("vector_completed_chunks")
              print(vector_completed_chunks)
              global_vector_chunks_to_be_spliced <<- vector_chunks_to_be_spliced
              print("vector_chunks_to_be_spliced")
              print(vector_chunks_to_be_spliced)
              
              global_list_result <<- list_result
              
              tibble_diffs <- tibble("n" = c(vector_current_chunks_spliced, vector_chunks_to_be_spliced) %>% sort, "n_minus_1" = c(c(vector_current_chunks_spliced, vector_chunks_to_be_spliced) %>% sort %>% .[2:length(.)], NA)) %>% 
                dplyr::mutate("diff" = `n_minus_1` - `n`)
              
              index_last_consecutive <- intersect(which(tibble_diffs$diff == 1), c(which(tibble_diffs$diff > 1) - 1, nrow(tibble_diffs) - 1)) %>% .[1]
              
              if (length(index_last_consecutive) > 0) {
                vector_chunks_to_be_spliced <- intersect(vector_chunks_to_be_spliced, tibble_diffs$n_minus_1 %>% .[1:index_last_consecutive])
              } else {
                vector_chunks_to_be_spliced <- integer(0)
              }
              
            }
            
          }
          
          # splice in terms
          if (length(vector_chunks_to_be_spliced) > 0) {
            
            list_new_chunks <- purrr::map(
              .x = vector_chunks_to_be_spliced,
              .f = function(a1) {
                
                chunk <- readRDS(file = paste(.temp_dir, "/", .job_name, "_chunk_", a1, ".rds", sep = ""))
                
                return(chunk)
                
              } )
            
            list_result <- purrr::splice(list_result, list_new_chunks)
            
            vector_current_chunks_spliced <- c(vector_current_chunks_spliced, vector_chunks_to_be_spliced) %>% sort
            
          }
          
          Sys.sleep(1)
          
        }
        
      } # END SPLICER ###
      
      # an actually useful progress bar although rudimentary
      if (ls(pattern = "^list_result$") %>% length == 0) {
        cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%)", sep = ""))
      } else {
        cat(paste("\r[", date(), "] Percent map completion: ", length(which(vector_logical_indices_workers_completed_reported)), "/", current_map_index, "/", map_length, " (", round(x = 100*length(which(vector_logical_indices_workers_completed_reported))/map_length, digits = 1), "%); Splice progress: ", length(list_result), "/", map_length, sep = ""))
      }
      
      # deal with completion flag
      if (.result_mode %in% c("ordered", "unordered")) {
        flag_completion <- length(list_result) == map_length
      } else {
        flag_completion <- length(which(vector_logical_indices_workers_completed_reported)) == map_length
      }
      
    }
    
    Sys.sleep(1)
    
  }
  
  # if any error, stop all
  if (any(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)) {
    purrr::map(.x = list_workers, .f = ~.x$signal(9))
    stop(print(paste("Exit status failure received on chunks:", paste(names(list_workers)[which(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)], collapse = ", "))))
  }
  
  # after the while loop is done, deal with splicing and socket transfer if we want the object returned immediately
  if (.result_mode %in% c("ordered", "unordered")) {
    return(list_result)
  } else {
    return(NULL)
  }
  
  if (.keep_intermediate_files == FALSE) {
    file.remove(list.files(path = .temp_dir, pattern = paste(.job_name, "_chunk_.*.rds", sep = ""), full.names = TRUE ))
  }
  
  rm(list = ls(pattern = "chunk_"))
  rm(list_workers)
  
  cat("\n")
  
}

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
  
  global_query_chr <<- query_chr
  global_query_start <<- query_start
  global_query_end <<- query_end
  global_query_strand <<- query_strand
  
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

extract_overlapping_features <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, left_query_shift = 0, right_query_shift = 0, left_tolerance = 0, right_tolerance = 0, return_type = NULL) {
  
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
split_delimited_columns_in_table <- function(input_table, target_colnames, split, columns_to_deduplicate = NULL) {
  
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
  
  list_split_lengths <- pmap(
      .l = list_split_columns, 
      .f = function(...) {
        return(unique(unlist(purrr::map(.x = list(...), .f = ~length(.x)))))
      }
    )
  
  if (unique(unlist(purrr::map(.x = list_split_lengths, .f = ~length(.x))))[1] != 1 | length(unique(unlist(purrr::map(.x = list_split_lengths, .f = ~length(.x))))) != 1) {
    stop("Multicolumn splits are uneven. Cannot proceed further.")
  }
  
  # repeat table according to the split lengths
  row_indices_of_table_repeated_by_split <- rep(x = 1:nrow(input_table), times = unlist(list_split_lengths))
  
  input_table_repeated_by_split <- input_table[row_indices_of_table_repeated_by_split, ]
  
  # replace target column with split values
  input_table_repeated_by_split[, target_colnames] <- input_table_repeated_by_split
  
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
    list_deduplicated_columns_split <- purrr::map(.x = list_deduplicated_columns, .f = function(a1) {
      
      # map a subset each of the L2 (elements of a column)
      a1[indices_of_duplicates] <- purrr::map2(.x = a1[indices_of_duplicates], .y = repetition_numbers_of_duplicates, 
                                               .f = ~rep(.x, times = .y) %>% unlist %>% paste(., 1:.y, sep = "_"))
      
      return(a1 %>% unlist)
      
    } )
    
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

# FUNCTION TO PLOT UMAP FOR SAMPLE AND REPLICATE
# tibble_metadata: this should have n columns equal to the number of annotation to be specified, and m rows equal to number of samples (including replicates). basically a collection of vectors of length equal to number of samples.
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
  
  transposed_matrixtable <- table_matrix %>% t
  # rownames(transposed_matrixtable) <- colnames(table_matrix)
  
  tibble_umap_result <- umap::umap(transposed_matrixtable) %>% .$layout %>% 
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

# FUNCTION TO PLOT tSNE FOR SAMPLE AND REPLICATE
# tibble_metadata: this should have n columns equal to the number of annotation to be specified, and m rows equal to number of samples (including replicates). basically a collection of vectors of length equal to number of samples.
# protected colname: "condition_names", "replicate_names"
plot_tSNE_for_timepoint_and_replicate_plotly <- function(table_matrix, timepoint_order = NULL, tibble_metadata = NULL, replicate_order = NULL, plot_shapes = TRUE, centroid_labels = TRUE, point_size = 5, centroid_label_size = 4, legend_position = "none", PCA_depths_y = NULL, PCA_depths_x = NULL, input_colour_limits = NULL, input_colour_value = NULL, save_dir = NULL, save_name = NULL, graph_title = NULL, width = 10, height = 10, ...) {
  
  # DEBUG ###
  # table_matrix <- tibble_matrix_absolute_psi_in_sample_replicate_format_with_na
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
  # tibble_metadata <- tibble_sharp_cluster_mapping %>% dplyr::ungroup() %>% dplyr::rename("som_cluster" = "cluster", "condition_names" = "condition_names") %>% dplyr::select(som_cluster, condition_names, Description, Biological_source, Type) %>% dplyr::mutate_at(.vars = c("Description", "Biological_source", "Type"), .funs = function(x) {gsub(x = x, pattern = "[^a-zA-Z0-9 -/]", replacement = "")}) 
  # input_colour_value = rainbow(tibble_sharp_cluster_mapping$cluster %>% max) %>% .[tibble_sharp_cluster_mapping$cluster]
  # save_dir <- R_processing_results_dir
  # save_name <- paste("atlas_polya_psisigma_consensus_som_clusterone_with_na_plotly", sep = "")
  # graph_title = "PolyA/PSI-Sigma w/ replicates, coloured by consensus SOM clusters"
  # width <- 40
  # height <- 40
  ###########
  
  transposed_matrixtable <- table_matrix %>% t
  
  tibble_tsne_result <- Rtsne::Rtsne(X = transposed_matrixtable %>% unique, perplexity = 1, check_duplicates = FALSE, verbose = TRUE, num_threads = 0) %>% .$Y %>% 
    as_tibble(.name_repair = "unique") %>%
    add_column("condition_names" = rownames(transposed_matrixtable %>% unique), .before = 1) %>%
    setNames(nm = c("condition|replicate", "V1", "V2")) %>%
    add_column("condition_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\1"),
               "replicate_names" = gsub(x = .$`condition|replicate`, pattern = "(.*)\\|(.*)", replacement = "\\2")
    ) %>% 
    dplyr::inner_join(., tibble_metadata, by = intersect(colnames(tibble_metadata), colnames(.))[intersect(colnames(tibble_metadata), colnames(.)) %in% c("condition_names", "replicate_names")])
  
  # centroid labels: make a separate file containing the centroid locations
  tibble_centroid_locations <- tibble_tsne_result %>% 
    dplyr::group_by(condition_names) %>% 
    dplyr::summarise("centroid_x" = mean(V1),
                     "centroid_y" = mean(V2)) %>%
    add_column("cluster_number" = paste(1:nrow(.)), .after = "condition_names")
  
  # append cluster info
  tibble_tsne_result_plot <- dplyr::left_join(tibble_tsne_result, tibble_centroid_locations, by = "condition_names")
  
  if (is.null(input_colour_limits) == TRUE) {
    
    input_colour_limits <- c(tibble_tsne_result_plot$condition_names %>% unique)
    
  }
  
  if (is.null(input_colour_value) == TRUE) {
    
    input_colour_value <- rainbow(n = (timepoint_order %>% length))
    
  }
  
  ggplot_plot <- ggplot() +
    (if (plot_shapes == TRUE) {geom_point(aes(x = tibble_tsne_result_plot$V1, y = tibble_tsne_result_plot$V2, shape = tibble_tsne_result_plot$replicate_names, color = tibble_tsne_result_plot$condition_names, fill = tibble_tsne_result_plot$condition_names), size = point_size) } else {geom_blank(aes(x = tibble_tsne_result_plot$V1, y = tibble_tsne_result_plot$V2))}) +
    scale_shape_manual(name = "Replicate", values = 1:length(replicate_order)) +
    (if (centroid_labels == TRUE) {geom_text(data = tibble_tsne_result_plot, mapping = do.call(what = aes, args = list("x" = tibble_tsne_result_plot$centroid_x, "y" = tibble_tsne_result_plot$centroid_y, "label" = tibble_tsne_result_plot$cluster_number, "color" = tibble_tsne_result_plot$condition_names) %>% purrr::splice(tibble_umap_result %>% dplyr::select(-`condition|replicate`, -V1, -V2, -contains("condition_names"), -contains("replicate_names")) %>% purrr::array_tree(margin = 2))), size = centroid_label_size)} else {geom_blank(aes(x = tibble_tsne_result_plot$centroid_x, y = tibble_tsne_result_plot$centroid_y))}) +
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

