source("/PHShome/kv1026/p05/remap_functions.R")

odir = paste0(imed, "ucsc_qtl_tracks/")
if (!dir.exists(odir)) dir.create(odir, recursive = TRUE)

load(paste0(imed, "Frank/qtls_list_260204.rda"))

qtls_name_map = c(
  "Brem_2002" = "Brem_2002",
  "Brem_2005" = "Brem_2005",
  "Smith_2008_Glucose" = "SmithGlc_2008",
  "Smith_2008_Ethanol" = "SmithEth_2008",
  "Khan_2009" = "Khan_2009",
  "Foss_2011" = "Foss_2011",
  "Picotti_2013" = "Picotti_2013",
  "Albert_2018" = "Albert_2018",
  "Grossbach_2022_mRNA" = "GrossbachRNA_2022",
  "Grossbach_2022_protein" = "GrossbachPRT_2022",
  "Albert_2014" = "Albert_2014",
  "Brion_2020_mRNA" = "BrionRNA_2020",
  "Brion_2020_protein" = "BrionPRT_2020"
)

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

format_chrom = function(x) {
  x = as.character(x)
  x = trimws(x)

  x = ifelse(grepl("^chr", x), x, paste0("chr", x))
  x
}

get_nice_name = function(study) {
  if (!(study %in% names(studyNiceNames))) {
    stop("Missing pretty name for study: ", study)
  }

  studyNiceNames[[study]]
}

rescale01 = function(x) {
  rng = range(x, na.rm = TRUE, finite = TRUE)

  if (!all(is.finite(rng))) {
    return(rep(1, length(x)))
  }

  if (diff(rng) == 0) {
    return(rep(1, length(x)))
  }

  (x - rng[1]) / diff(rng)
}

blend_color = function(low, high, t) {
  t = pmin(pmax(t, 0), 1)
  low_rgb = grDevices::col2rgb(low)
  high_rgb = grDevices::col2rgb(high)

  vapply(
    t,
    function(z) {
      rgb = round((1 - z) * low_rgb + z * high_rgb)
      paste(rgb[1, 1], rgb[2, 1], rgb[3, 1], sep = ",")
    },
    character(1)
  )
}

make_track = function(name, description) {
  paste(
    paste0("track name=\"", name, "\""),
    paste0("description=\"", description, "\""),
    "type=bed",
    "visibility=dense",
    "itemRgb=On",
    sep = " "
  )
}

required_cols = c(
  "gene_ID",
  "marker_name",
  "cor_coef",
  "LOD",
  "CI_min_base",
  "CI_max_base",
  "CI_chr"
)

missing_lists = setdiff(names(qtls_name_map), names(qtls_list))
if (length(missing_lists) > 0L) {
  stop("Missing expected qtls_list entries: ", paste(missing_lists, collapse = ", "))
}

qtls = lapply(names(qtls_name_map), function(list_name) {
  df = qtls_list[[list_name]]
  missing_cols = setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns for ",
      list_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df %>%
    dplyr::mutate(
      source_name = list_name,
      study = unname(qtls_name_map[[list_name]])
    )
}) %>%
  dplyr::bind_rows()

qtls = qtls %>%
  dplyr::filter(study != "AlbertAlt_2018")

peak = parse_marker(qtls$marker_name)

qtls$chrom = format_chrom(qtls$CI_chr)

qtls$start = pmin(qtls$CI_min_base, qtls$CI_max_base, na.rm = TRUE)
qtls$end = pmax(qtls$CI_min_base, qtls$CI_max_base, na.rm = TRUE)

missing_interval = !is.finite(qtls$start) | !is.finite(qtls$end)
qtls$start[missing_interval] = peak$pos[missing_interval]
qtls$end[missing_interval] = peak$pos[missing_interval]

qtls$chromStart = pmax(as.integer(qtls$start) - 1L, 0L)
qtls$chromEnd = as.integer(qtls$end)

qtls$thickStart = pmax(as.integer(peak$pos) - 1L, qtls$chromStart)
qtls$thickEnd = pmin(as.integer(peak$pos), qtls$chromEnd)

bad_width = !is.finite(qtls$chromEnd) | qtls$chromEnd <= qtls$chromStart
qtls$chromEnd[bad_width] = qtls$chromStart[bad_width] + 1L

bad_thick = !is.finite(qtls$thickEnd) | qtls$thickEnd <= qtls$thickStart
qtls$thickStart[bad_thick] = qtls$chromStart[bad_thick]
qtls$thickEnd[bad_thick] = qtls$chromStart[bad_thick] + 1L

qtls$thickEnd = pmin(qtls$thickEnd, qtls$chromEnd)

qtls$strand = ifelse(qtls$cor_coef < 0, "-", "+")

qtls = qtls %>%
  dplyr::group_by(study) %>%
  dplyr::mutate(
    lod_scaled = rescale01(LOD),
    itemRgb = ifelse(
      cor_coef < 0,
      blend_color(low = "#DEEBF7", high = "#08519C", t = lod_scaled),
      blend_color(low = "#FEE0D2", high = "#DE2D26", t = lod_scaled)
    ),
    score = pmin(1000L, pmax(0L, as.integer(round(1000 * lod_scaled))))
  ) %>%
  dplyr::ungroup()

qtls$name = paste0(
  qtls$gene_ID,
  "_",
  qtls$marker_name,
  "_LOD",
  sprintf("%.2f", qtls$LOD)
)

bed = qtls[, c(
  "chrom", "chromStart", "chromEnd", "name", "score",
  "strand", "thickStart", "thickEnd", "itemRgb"
)]

outfile = paste0(odir, "frank_qtls_ucsc_multitrack.bed")
con = file(outfile, "wt")

writeLines("browser db=sacCer3", con)

study_order = unname(qtls_name_map)
study_order = study_order[study_order %in% unique(qtls$study)]

for (st in study_order) {
  nice = get_nice_name(st)

  writeLines(
    make_track(
      name = nice,
      description = paste0(nice, " QTLs")
    ),
    con
  )

  df = bed[qtls$study == st, , drop = FALSE]
  df = df[order(df$chrom, df$chromStart, df$chromEnd, df$name), , drop = FALSE]

  write.table(
    df,
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

close(con)

message("Wrote UCSC multitrack BED to: ", outfile)