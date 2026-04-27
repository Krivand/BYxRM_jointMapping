source("/PHShome/kv1026/p05/remap_functions.R")

load(paste0(imed, "indivFullQTLs.rda"))

parse_marker = function(x) {
  x = as.character(x)

  chrom = sub(":.*$", "", x)
  pos = suppressWarnings(as.integer(sub("^.*:([0-9]+).*$", "\\1", x)))

  data.frame(
    chrom = chrom,
    pos = pos,
    stringsAsFactors = FALSE
  )
}

cat("========================================\n")
cat("TRACK DIAGNOSTICS\n")
cat("========================================\n\n")

cil  = parse_marker(indivFullQTLs$CI.l)
cir  = parse_marker(indivFullQTLs$CI.r)
peak = parse_marker(indivFullQTLs$pmarker)

# ---- NA diagnostics ----
cat("---- NA COUNTS ----\n")
cat("CI.l  NA pos:", sum(is.na(cil$pos)), "\n")
cat("CI.r  NA pos:", sum(is.na(cir$pos)), "\n")
cat("Peak  NA pos:", sum(is.na(peak$pos)), "\n\n")

# ---- Examples of bad rows ----
show_examples = function(vec, mask, label) {
  if (sum(mask) == 0) {
    cat(label, ": none\n")
    return()
  }
  cat(label, " (showing up to 10):\n")
  print(unique(vec[mask])[1:min(10, sum(mask))])
  cat("\n")
}

show_examples(indivFullQTLs$CI.l, is.na(cil$pos),  "Bad CI.l")
show_examples(indivFullQTLs$CI.r, is.na(cir$pos),  "Bad CI.r")
show_examples(indivFullQTLs$pmarker, is.na(peak$pos), "Bad pmarker")

# ---- Cross chromosome intervals ----
cross_chr = !is.na(cil$chrom) & !is.na(cir$chrom) & cil$chrom != cir$chrom

cat("---- CROSS-CHROMOSOME CIs ----\n")
cat("Count:", sum(cross_chr), "\n")

if (sum(cross_chr) > 0) {
  print(head(indivFullQTLs[cross_chr, c("study","gene","CI.l","CI.r")], 10))
}
cat("\n")

# ---- Coordinate sanity ----
start = pmin(cil$pos, cir$pos, na.rm = TRUE)
end   = pmax(cil$pos, cir$pos, na.rm = TRUE)

bad_order = is.finite(start) & is.finite(end) & (start > end)

cat("---- START > END ISSUES ----\n")
cat("Count:", sum(bad_order), "\n")

if (sum(bad_order) > 0) {
  print(head(indivFullQTLs[bad_order, c("study","gene","CI.l","CI.r")], 10))
}
cat("\n")

# ---- Zero-width intervals ----
zero_width = is.finite(start) & is.finite(end) & (start == end)

cat("---- ZERO-WIDTH (1bp) INTERVALS ----\n")
cat("Count:", sum(zero_width), "\n\n")

# ---- Missing chrom labels ----
missing_chr = is.na(cil$chrom) | is.na(cir$chrom)

cat("---- MISSING CHROMOSOMES ----\n")
cat("Count:", sum(missing_chr), "\n")

if (sum(missing_chr) > 0) {
  print(head(indivFullQTLs[missing_chr, c("study","gene","CI.l","CI.r")], 10))
}
cat("\n")

cat("========================================\n")
cat("DONE\n")
cat("========================================\n")
