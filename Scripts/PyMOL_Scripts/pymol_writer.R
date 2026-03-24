write_pymol_colouring_script <- function() {
  
  # ---- pull required globals ----
  cutoff       <- get("cutoff",       envir = .GlobalEnv)
  cutoffstrict <- get("cutoffstrict", envir = .GlobalEnv)
  protein      <- get("protein",         envir = .GlobalEnv)
  PREDNA       <- get("PREDNA",       envir = .GlobalEnv)
  
  # ---- resolve paths using here() ----
  data2bfactor_path <- here::here(
    "Scripts", "PyMOL_Scripts", "data2bfactor.py"
  )
  
  spectrumany_path <- here::here(
    "Scripts", "PyMOL_Scripts", "spectrumany.py"
  )
  
  ratio_file_path <- here::here(
    "Output", "PyMOL", "Ratio_Files",
    paste0(protein, "_", PREDNA, "_scaled.txt")
  )
  
  # ---- output location ----
  out_dir <- here::here("PyMOL_Structures", protein)
  out_file <- here::here(
    "PyMOL_Structures",
    protein,
    paste0(protein, "_", PREDNA, "_colouring.txt")
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- build PyMOL script (ALL colours static) ----
  script_text <- sprintf("
run %s
run %s

set use_shaders, 1

alter prot, b = -1
data2b_res prot, %s

color grey50, prot

spectrumany b, 0x7c8fb3 0x858585 0x737373, prot, %f, 1

select broad, prot and ((b > %f and b < %f) or b = %f)
spectrumany b, 0x65b4bf 0x648bcc 0x7c8fb3, broad, %f, %f

select strict, prot and (b > -1 and b < %f)
spectrumany b, 0x9eda4b 0x7ec46a 0x65b4bf, strict, 0, %f

select unassigned, prot and b = -1
color 0xd9d9d9, unassigned

set ray_trace_mode, 1
set ambient, 0.5
set cartoon_nucleic_acid_color, 0x8046b2
set cartoon_ring_mode, 0
set cartoon_ladder_color, 0xb274de
set cartoon_ring_color, 0xb274de
set bg_rgb, white
set ray_shadows, 1
",
                         data2bfactor_path,
                         spectrumany_path,
                         ratio_file_path,
                         cutoff,
                         cutoffstrict, cutoff, cutoff,
                         cutoffstrict, cutoff,
                         cutoffstrict,
                         cutoffstrict
  )
  
  writeLines(script_text, out_file)
  
  invisible(out_file)
}