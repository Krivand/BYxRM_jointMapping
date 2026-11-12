library(tidyverse)

inps = c(
  "/PHShome/kv1026/p05/imed/jointAssocsFDR_geno",
  "/PHShome/kv1026/p05/imed/jointAssocsFDR_analyte"
)

outs = c(
  "/PHShome/kv1026/p05/imed/jointAnalyte_geno.rda",
  "/PHShome/kv1026/p05/imed/jointAnalyte_ixn.rda"
)

stopifnot(length(inps) == length(outs))

aggregate_one = function(inp_dir, out_rda) {
  files = list.files(inp_dir, pattern = "^chunk[0-9]+\\.rda$", full.names = TRUE)
  if (length(files) == 0) stop("No chunk*.rda in: ", inp_dir)

  # numeric sort by chunk index
  idx = as.integer(sub("^chunk([0-9]+)\\.rda$", "\\1", basename(files)))
  files = files[order(idx)]

  agg = NULL
  for (f in files) {
    e = new.env(parent = emptyenv())
    load(f, envir = e)
    x = e$chunkAssocs
    if (is.null(x)) stop("chunkAssocs missing in: ", f)

    if (is.null(agg)) agg = x
    else agg = dplyr::bind_rows(agg, x)
  }

  out_name = tools::file_path_sans_ext(basename(out_rda))
  assign(out_name, agg)

  out_dir = dirname(out_rda)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  save(list = out_name, file = out_rda)
  message("Wrote ", out_rda, " (object: ", out_name, ", rows: ", nrow(agg), ")")
}

for (i in seq_along(inps)) {
  aggregate_one(inps[[i]], outs[[i]])
}
