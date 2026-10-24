
cat("============================================\n")
cat("SCRIPT 07b: TRACE-N-SEQ METADATA AUDIT\n")
cat("============================================\n\n")

outdir <- "data/raw/E-MTAB-12940"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

base_url <- paste0(
  "https://ftp.ebi.ac.uk/biostudies/fire/",
  "E-MTAB-/940/E-MTAB-12940/Files/"
)

files <- c(
  "smartseq_metadata_table.txt",
  "smartseq_metadata_table2.txt"
)

## ---------------------------------------------------------
## Download metadata files
## ---------------------------------------------------------

for (f in files) {

  dest <- file.path(outdir, f)

  if (!file.exists(dest) || file.info(dest)$size == 0) {

    cat("Downloading:", f, "\n")

    download.file(
      paste0(base_url, f),
      destfile = dest,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    )

  } else {

    cat("Already present:", f, "\n")
  }
}

cat("\n============================================\n")
cat("FILE SIZES / MD5\n")
cat("============================================\n")

for (f in files) {

  path <- file.path(outdir, f)

  cat(
    f,
    " | bytes =", file.info(path)$size,
    " | MD5 =", unname(tools::md5sum(path)),
    "\n"
  )
}

## ---------------------------------------------------------
## Robust table reader
## ---------------------------------------------------------

read_meta <- function(path) {

  x <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )

  if (is.null(x)) {
    stop("Could not read: ", path)
  }

  x
}

m1 <- read_meta(file.path(outdir, files[1]))
m2 <- read_meta(file.path(outdir, files[2]))

## ---------------------------------------------------------
## Audit function
## ---------------------------------------------------------

audit_meta <- function(x, label) {

  cat("\n\n============================================\n")
  cat(label, "\n")
  cat("============================================\n")

  cat("Rows:", nrow(x), "\n")
  cat("Columns:", ncol(x), "\n")

  cat("\nCOLUMN NAMES\n")
  print(names(x))

  ## Find biologically relevant columns
  pattern <- paste(
    "cell",
    "sample",
    "disease",
    "condition",
    "organ",
    "tissue",
    "gang",
    "origin",
    "type",
    "annotation",
    "quality",
    "qc",
    "ident",
    "individual",
    "mouse",
    "sex",
    sep = "|"
  )

  candidate_cols <- names(x)[
    grepl(pattern, names(x), ignore.case = TRUE)
  ]

  cat("\nCANDIDATE METADATA COLUMNS\n")
  print(candidate_cols)

  ## Print manageable value distributions
  for (nm in candidate_cols) {

    vals <- as.character(x[[nm]])
    vals[is.na(vals)] <- "<NA>"

    tab <- sort(table(vals), decreasing = TRUE)

    cat("\n---", nm, "---\n")
    cat("Unique values:", length(tab), "\n")

    print(head(tab, 30))
  }

  ## Search rows containing terms central to our analysis
  terms <- c(
    "pancreas",
    "healthy",
    "DRG",
    "JNG",
    "CG",
    "sensory",
    "peptidergic",
    "neurofilament"
  )

  text_matrix <- apply(
    x,
    1,
    function(z) paste(z, collapse = " | ")
  )

  cat("\nKEYWORD ROW COUNTS\n")

  for (term in terms) {

    n <- sum(
      grepl(
        term,
        text_matrix,
        ignore.case = TRUE
      )
    )

    cat(term, ":", n, "\n")
  }

  invisible(candidate_cols)
}

cols1 <- audit_meta(
  m1,
  "smartseq_metadata_table.txt"
)

cols2 <- audit_meta(
  m2,
  "smartseq_metadata_table2.txt"
)

## ---------------------------------------------------------
## Compare versions
## ---------------------------------------------------------

cat("\n\n============================================\n")
cat("VERSION COMPARISON\n")
cat("============================================\n")

cat(
  "Rows old/new:",
  nrow(m1), "/",
  nrow(m2), "\n"
)

cat(
  "Columns old/new:",
  ncol(m1), "/",
  ncol(m2), "\n"
)

common_cols <- intersect(names(m1), names(m2))

cat(
  "Common columns:",
  length(common_cols), "\n"
)

cat(
  "Only metadata_table:",
  paste(setdiff(names(m1), names(m2)), collapse = "; "),
  "\n"
)

cat(
  "Only metadata_table2:",
  paste(setdiff(names(m2), names(m1)), collapse = "; "),
  "\n"
)

## ---------------------------------------------------------
## Search specifically for likely pancreas-projecting cells
## ---------------------------------------------------------

extract_keyword_rows <- function(x) {

  row_text <- apply(
    x,
    1,
    function(z) paste(z, collapse = " | ")
  )

  pancreas <- grepl(
    "pancreas",
    row_text,
    ignore.case = TRUE
  )

  drg <- grepl(
    "\\bDRG\\b",
    row_text,
    ignore.case = TRUE
  )

  jng <- grepl(
    "\\bJNG\\b",
    row_text,
    ignore.case = TRUE
  )

  sensory <- grepl(
    "sensory|peptidergic|neurofilament",
    row_text,
    ignore.case = TRUE
  )

  data.frame(
    pancreas = pancreas,
    DRG = drg,
    JNG = jng,
    sensory = sensory
  )
}

k1 <- extract_keyword_rows(m1)
k2 <- extract_keyword_rows(m2)

cat("\n============================================\n")
cat("PANCREAS / GANGLION KEYWORD INTERSECTIONS\n")
cat("============================================\n")

cat("\nmetadata_table.txt\n")
cat("Pancreas rows:", sum(k1$pancreas), "\n")
cat("Pancreas + DRG:", sum(k1$pancreas & k1$DRG), "\n")
cat("Pancreas + JNG:", sum(k1$pancreas & k1$JNG), "\n")
cat(
  "Pancreas + sensory annotation:",
  sum(k1$pancreas & k1$sensory),
  "\n"
)

cat("\nmetadata_table2.txt\n")
cat("Pancreas rows:", sum(k2$pancreas), "\n")
cat("Pancreas + DRG:", sum(k2$pancreas & k2$DRG), "\n")
cat("Pancreas + JNG:", sum(k2$pancreas & k2$JNG), "\n")
cat(
  "Pancreas + sensory annotation:",
  sum(k2$pancreas & k2$sensory),
  "\n"
)

## Save audit objects
saveRDS(
  list(
    metadata_table = m1,
    metadata_table2 = m2
  ),
  file.path(
    outdir,
    "E-MTAB-12940_metadata_audit.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  "reproducibility/07b_EMTAB12940_metadata_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 07b COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

