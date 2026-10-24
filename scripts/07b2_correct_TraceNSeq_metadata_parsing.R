
cat("============================================\n")
cat("SCRIPT 07b2: CORRECT TRACE-N-SEQ METADATA PARSING\n")
cat("============================================\n\n")

indir <- "data/raw/E-MTAB-12940"

f1 <- file.path(indir, "smartseq_metadata_table.txt")
f2 <- file.path(indir, "smartseq_metadata_table2.txt")

stopifnot(file.exists(f1), file.exists(f2))

## ============================================================
## READ THE TWO FILES USING THEIR ACTUAL DELIMITERS
## ============================================================

m1 <- read.delim(
    f1,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
)

m2 <- read.csv(
    f2,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
)

cat("Batch 1 dimensions:", nrow(m1), "x", ncol(m1), "\n")
cat("Batch 2 dimensions:", nrow(m2), "x", ncol(m2), "\n")

cat("\nBatch 1 columns:\n")
print(names(m1))

cat("\nBatch 2 columns:\n")
print(names(m2))

## ============================================================
## AUDIT FUNCTION
## ============================================================

audit_batch <- function(x, batch) {

    cat("\n\n============================================\n")
    cat(batch, "\n")
    cat("============================================\n")

    fields <- c(
        "Ganglion",
        "MouseModel",
        "FastBlue",
        "Tumor",
        "CellLine",
        "InnervatingTissue",
        "Replicate",
        "Is_Cell",
        "is_cell",
        "ClusterNew",
        "AnnotatedCellType",
        "AnnotatedSubCellType",
        "AnnotatedCellTypeSharma"
    )

    fields <- intersect(fields, names(x))

    for (nm in fields) {

        cat("\n---", nm, "---\n")

        tab <- sort(
            table(
                x[[nm]],
                useNA = "ifany"
            ),
            decreasing = TRUE
        )

        print(tab)
    }

    ## --------------------------------------------------------
    ## Ganglion × tissue
    ## --------------------------------------------------------

    if (all(c("Ganglion", "InnervatingTissue") %in% names(x))) {

        cat("\n===== GANGLION x INNERVATING TISSUE =====\n")

        print(
            table(
                Ganglion = x$Ganglion,
                Tissue = x$InnervatingTissue,
                useNA = "ifany"
            )
        )
    }

    ## --------------------------------------------------------
    ## Ganglion × mouse model
    ## --------------------------------------------------------

    if (all(c("Ganglion", "MouseModel") %in% names(x))) {

        cat("\n===== GANGLION x MOUSE MODEL =====\n")

        print(
            table(
                Ganglion = x$Ganglion,
                Model = x$MouseModel,
                useNA = "ifany"
            )
        )
    }

    ## --------------------------------------------------------
    ## Pancreas-labelled rows
    ## --------------------------------------------------------

    if ("InnervatingTissue" %in% names(x)) {

        pancreas <- tolower(trimws(x$InnervatingTissue)) == "pancreas"

        cat("\n===== PANCREAS-LABELLED ROWS =====\n")
        cat("Total:", sum(pancreas, na.rm = TRUE), "\n")

        if ("Ganglion" %in% names(x)) {
            cat("\nGanglion:\n")
            print(table(x$Ganglion[pancreas], useNA = "ifany"))
        }

        if ("MouseModel" %in% names(x)) {
            cat("\nMouse model:\n")
            print(table(x$MouseModel[pancreas], useNA = "ifany"))
        }

        if ("Tumor" %in% names(x)) {
            cat("\nTumor field:\n")
            print(table(x$Tumor[pancreas], useNA = "ifany"))
        }

        if ("CellLine" %in% names(x)) {
            cat("\nCellLine:\n")
            print(table(x$CellLine[pancreas], useNA = "ifany"))
        }

        if ("FastBlue" %in% names(x)) {
            cat("\nFastBlue:\n")
            print(table(x$FastBlue[pancreas], useNA = "ifany"))
        }

        if ("AnnotatedCellType" %in% names(x)) {
            cat("\nAnnotated cell type:\n")

            print(
                sort(
                    table(
                        x$AnnotatedCellType[pancreas],
                        useNA = "ifany"
                    ),
                    decreasing = TRUE
                )
            )
        }

        if ("AnnotatedSubCellType" %in% names(x)) {
            cat("\nAnnotated subtype:\n")

            print(
                sort(
                    table(
                        x$AnnotatedSubCellType[pancreas],
                        useNA = "ifany"
                    ),
                    decreasing = TRUE
                )
            )
        }
    }

    ## --------------------------------------------------------
    ## Search entire metadata for vagal/JNG terminology
    ## --------------------------------------------------------

    charcols <- names(x)[vapply(
        x,
        function(z) is.character(z) || is.factor(z),
        logical(1)
    )]

    if (length(charcols) > 0) {

        txt <- apply(
            x[, charcols, drop = FALSE],
            1,
            paste,
            collapse = " | "
        )

        terms <- c(
            "JNG",
            "jugular",
            "nodose",
            "vagus",
            "vagal"
        )

        cat("\n===== JNG / VAGAL TERMINOLOGY SEARCH =====\n")

        for (term in terms) {

            hit <- grepl(
                term,
                txt,
                ignore.case = TRUE
            )

            cat(term, ":", sum(hit), "\n")
        }
    }

    invisible(NULL)
}

audit_batch(m1, "SMARTSEQ BATCH 1")
audit_batch(m2, "SMARTSEQ BATCH 2")

## ============================================================
## COMBINED AUDIT
## ============================================================

cat("\n\n============================================\n")
cat("COMBINED BASIC COUNTS\n")
cat("============================================\n")

cat(
    "Total rows across both batches:",
    nrow(m1) + nrow(m2),
    "\n"
)

for (nm in c("Ganglion", "InnervatingTissue")) {

    if (nm %in% names(m1) && nm %in% names(m2)) {

        cat("\nCombined", nm, "counts:\n")

        vals <- c(
            as.character(m1[[nm]]),
            as.character(m2[[nm]])
        )

        print(
            sort(
                table(vals, useNA = "ifany"),
                decreasing = TRUE
            )
        )
    }
}

## ============================================================
## SAVE CORRECTLY PARSED OBJECTS
## ============================================================

saveRDS(
    m1,
    file.path(
        indir,
        "smartseq_metadata_batch1_parsed.rds"
    )
)

saveRDS(
    m2,
    file.path(
        indir,
        "smartseq_metadata_batch2_parsed.rds"
    )
)

saveRDS(
    list(
        batch1 = m1,
        batch2 = m2
    ),
    file.path(
        indir,
        "E-MTAB-12940_metadata_correctly_parsed.rds"
    )
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07b2_EMTAB12940_metadata_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 07b2 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

