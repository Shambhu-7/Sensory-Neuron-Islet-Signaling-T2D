
options(stringsAsFactors = FALSE)
options(timeout = 600)

outdir <- file.path(
    getwd(),
    "results",
    "ligands"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

base <- "https://omnipathdb.org"

# ============================================================
# Robust downloader
# ============================================================

download_table <- function(endpoint, query, outfile,
                           attempts = 3) {

    u <- paste0(
        base,
        endpoint,
        "?",
        query
    )

    tmp <- tempfile(fileext = ".tsv")

    for (attempt in seq_len(attempts)) {

        cat(
            "\nAttempt",
            attempt,
            "of",
            attempts,
            "\n"
        )

        status <- try(
            download.file(
                u,
                destfile = tmp,
                method = "libcurl",
                mode = "wb",
                quiet = FALSE
            ),
            silent = TRUE
        )

        if (
            !inherits(status, "try-error") &&
            file.exists(tmp) &&
            file.size(tmp) > 0
        ) {
            break
        }

        if (attempt == attempts) {
            stop(
                "Download failed after ",
                attempts,
                " attempts:\n",
                u
            )
        }

        Sys.sleep(5)
    }

    x <- read.delim(
        tmp,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    write.table(
        x,
        outfile,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat(
        "Rows:",
        nrow(x),
        "| Columns:",
        ncol(x),
        "\n"
    )

    x
}

# ============================================================
# Curated cell-cell communication resources
# ============================================================

curated_resources <- c(
    "CellPhoneDB",
    "Guide2Pharma",
    "HPMR",
    "ICELLNET",
    "Kirouac2010",
    "CellTalkDB",
    "CellChatDB",
    "connectomeDB2020",
    "Ramilowski2015",
    "talklr"
)

# ============================================================
# Download each interaction resource separately
#
# This is deliberately preferable to downloading the complete
# omnipath + ligrecextra network.
# ============================================================

interaction_list <- list()
resource_status <- list()

for (res in curated_resources) {

    cat("\n============================================\n")
    cat("RESOURCE:", res, "\n")
    cat("============================================\n")

    query <- paste0(
        "genesymbols=yes",
        "&datasets=omnipath,ligrecextra",
        "&organisms=9606",
        "&databases=", res,
        "&fields=sources,references,curation_effort",
        "&license=academic"
    )

    outfile <- file.path(
        outdir,
        paste0(
            "OmniPath_DIRECT_interactions_",
            res,
            ".tsv"
        )
    )

    ans <- try(
        download_table(
            "/interactions",
            query,
            outfile
        ),
        silent = TRUE
    )

    if (inherits(ans, "try-error")) {

        cat(
            "WARNING: resource retrieval failed:",
            res,
            "\n"
        )

        resource_status[[res]] <- data.frame(
            resource = res,
            status = "FAILED",
            rows = NA_integer_,
            stringsAsFactors = FALSE
        )

        next
    }

    if (nrow(ans) == 0) {

        cat(
            "WARNING: zero records:",
            res,
            "\n"
        )

        resource_status[[res]] <- data.frame(
            resource = res,
            status = "ZERO_ROWS",
            rows = 0L,
            stringsAsFactors = FALSE
        )

        next
    }

    ans$query_resource <- res

    interaction_list[[res]] <- ans

    resource_status[[res]] <- data.frame(
        resource = res,
        status = "SUCCESS",
        rows = nrow(ans),
        stringsAsFactors = FALSE
    )
}

# ============================================================
# Resource retrieval status
# ============================================================

resource_status <- do.call(
    rbind,
    resource_status
)

rownames(resource_status) <- NULL

write.csv(
    resource_status,
    file.path(
        outdir,
        "OmniPath_DIRECT_resource_status.csv"
    ),
    row.names = FALSE
)

cat("\n============================================\n")
cat("RESOURCE RETRIEVAL STATUS\n")
cat("============================================\n")

print(
    resource_status,
    row.names = FALSE
)

# ============================================================
# Combine successful interaction tables
# ============================================================

if (length(interaction_list) == 0) {
    stop("No interaction resources were retrieved.")
}

all_names <- unique(
    unlist(
        lapply(
            interaction_list,
            names
        )
    )
)

interaction_list <- lapply(
    interaction_list,
    function(x) {

        missing_cols <- setdiff(
            all_names,
            names(x)
        )

        for (z in missing_cols) {
            x[[z]] <- NA
        }

        x[
            ,
            all_names,
            drop = FALSE
        ]
    }
)

interactions <- do.call(
    rbind,
    interaction_list
)

rownames(interactions) <- NULL

write.table(
    interactions,
    file.path(
        outdir,
        "OmniPath_DIRECT_curated_interactions_combined.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\nCombined interaction records:",
    nrow(interactions),
    "\n"
)

# ============================================================
# Validate essential interaction columns
# ============================================================

needed_int <- c(
    "source",
    "target",
    "source_genesymbol",
    "target_genesymbol"
)

cat("\n============================================\n")
cat("INTERACTION ESSENTIAL COLUMNS\n")
cat("============================================\n")

essential_check <- data.frame(
    column = needed_int,
    present = needed_int %in% names(interactions)
)

print(
    essential_check,
    row.names = FALSE
)

if (!all(essential_check$present)) {

    cat(
        "\nActual interaction columns:\n"
    )

    print(
        names(interactions)
    )

    stop(
        "Required interaction columns missing."
    )
}

# ============================================================
# Ligand annotations
#
# Use generic/composite functional ligand classification.
# ============================================================

cat("\n============================================\n")
cat("DOWNLOADING LIGAND ANNOTATIONS\n")
cat("============================================\n")

ligand_query <- paste0(
    "categories=ligand",
    "&scope=generic",
    "&aspect=functional",
    "&source=composite",
    "&causality=trans",
    "&license=academic"
)

ligands <- download_table(
    "/intercell",
    ligand_query,
    file.path(
        outdir,
        "OmniPath_DIRECT_ligand_annotations_RAW.tsv"
    )
)

# ============================================================
# Receptor annotations
# ============================================================

cat("\n============================================\n")
cat("DOWNLOADING RECEPTOR ANNOTATIONS\n")
cat("============================================\n")

receptor_query <- paste0(
    "categories=receptor",
    "&scope=generic",
    "&aspect=functional",
    "&source=composite",
    "&causality=rec",
    "&license=academic"
)

receptors <- download_table(
    "/intercell",
    receptor_query,
    file.path(
        outdir,
        "OmniPath_DIRECT_receptor_annotations_RAW.tsv"
    )
)

# ============================================================
# Inspect intercell columns
# ============================================================

cat("\n============================================\n")
cat("LIGAND ANNOTATION COLUMNS\n")
cat("============================================\n")

print(
    names(ligands)
)

cat("\nLigand rows:", nrow(ligands), "\n")

cat("\n============================================\n")
cat("RECEPTOR ANNOTATION COLUMNS\n")
cat("============================================\n")

print(
    names(receptors)
)

cat("\nReceptor rows:", nrow(receptors), "\n")

# ============================================================
# Preview
# ============================================================

cat("\n============================================\n")
cat("LIGAND PREVIEW\n")
cat("============================================\n")

print(
    head(
        ligands,
        5
    )
)

cat("\n============================================\n")
cat("RECEPTOR PREVIEW\n")
cat("============================================\n")

print(
    head(
        receptors,
        5
    )
)

# ============================================================
# Save retrieval metadata
# ============================================================

metadata <- data.frame(
    retrieval_time = as.character(Sys.time()),
    timeout_seconds = getOption("timeout"),
    interaction_records = nrow(interactions),
    ligand_annotation_records = nrow(ligands),
    receptor_annotation_records = nrow(receptors),
    successful_resources =
        sum(resource_status$status == "SUCCESS"),
    failed_resources =
        sum(resource_status$status == "FAILED"),
    zero_row_resources =
        sum(resource_status$status == "ZERO_ROWS"),
    stringsAsFactors = FALSE
)

write.csv(
    metadata,
    file.path(
        outdir,
        "OmniPath_DIRECT_retrieval_metadata.csv"
    ),
    row.names = FALSE
)

cat("\n============================================\n")
cat("DIRECT OMNIPATH API TEST COMPLETED\n")
cat("============================================\n")

