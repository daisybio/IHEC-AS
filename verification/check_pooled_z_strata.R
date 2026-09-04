## Does pooling z ACROSS Event Types actually hold up for the minority type?
##
## SE is 95% of events, so a fully pooled reference is ~an SE reference, and an
## aggregate calibration number is dominated by SE. Test calibration SEPARATELY
## for SE and RI values, under three reference choices:
##   (a) reference pooled over ALL events (both types)
##   (b) reference pooled within Event Type   <- matches the FDR family
##   (c) reference pooled within Event Type x n-tertile
## Also report the resolution each choice would have AT FULL SCALE, since the
## earlier stratified test was resolution-limited only by the sample size.
setwd("/nfs/proj/quirinmanz/ihec_as_clean/IHEC-AS")
suppressPackageStartupMessages(library(data.table))

D <- "processed_data/event_models/biotype_filtered/screen"
set.seed(42)
# RI has few events, so take as many as possible for it
N_SE <- 700L
N_RI <- 700L

scr <- fread("processed_data/event_models/biotype_filtered/screen_results.csv.gz")[
  is.finite(p_emp) & R_used >= 150L
]
scr[, ID := as.integer(ID)]
sess <- readRDS("processed_data/session_09_1_ml_local_biotype_filtered.rds")
ed <- sess$event_dt[, .(ID = as.integer(ID), ET = `Event Type`)]
rm(sess); gc()
scr[ed, on = "ID", ET := i.ET]
scr[, n_tert := cut(n_samples, quantile(n_samples, 0:3 / 3), include.lowest = TRUE,
                    labels = c("small_n", "mid_n", "large_n")), by = ET]
samp <- rbind(
  scr[ET == "SE"][sample(.N, min(.N, N_SE))],
  scr[ET == "RI"][sample(.N, min(.N, N_RI))]
)
cat(sprintf("sampled SE=%d RI=%d\n", sum(samp$ET == "SE"), sum(samp$ET == "RI")))

read_z <- function(id) {
  f <- file.path(D, sprintf("%d_screen_null.csv.gz", id))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  if (!"feature_set" %in% names(d)) d[, feature_set := "long"]
  d <- d[feature_set == "long" & is.finite(null_R2)]
  if (nrow(d) < 100L) return(NULL)
  x <- d$null_R2; n <- length(x)
  s <- sum(x); ss <- sum(x^2)
  m_i <- (s - x) / (n - 1)
  v_i <- (ss - x^2 - (n - 1) * m_i^2) / (n - 2)
  data.table(ID = id, z = (x - m_i) / sqrt(pmax(v_i, .Machine$double.eps)))
}
.p <- pbmcapply::pbmclapply(samp$ID, read_z, mc.cores = 8)
nl <- rbindlist(.p[vapply(.p, is.data.frame, logical(1))], fill = TRUE)
nl[samp, on = "ID", `:=`(ET = i.ET, n_tert = i.n_tert)]
cat(sprintf("null z values: SE=%d RI=%d\n\n",
            nrow(nl[ET == "SE"]), nrow(nl[ET == "RI"])))

# leave-one-EVENT-out p against a given reference set of z
loo_p <- function(z_query, id_query, z_ref, id_ref) {
  R <- sort(z_ref); N <- length(R)
  ge <- N - findInterval(z_query - 1e-12, R)
  own <- integer(length(z_query))
  n_own <- integer(length(z_query))
  for (g in unique(id_query)) {
    iq <- which(id_query == g)
    zr <- z_ref[id_ref == g]
    if (length(zr)) {
      zz <- sort(zr)
      own[iq] <- length(zr) - findInterval(z_query[iq] - 1e-12, zz)
      n_own[iq] <- length(zr)
    }
  }
  (1 + (ge - own)) / (1 + (N - n_own))
}

ALPHAS <- c(1e-2, 6.64e-4, 2.8e-4, 1e-4)
show <- function(p, label) {
  p <- p[is.finite(p)]
  cat(sprintf("\n%-46s n=%d\n", label, length(p)))
  print(rbindlist(lapply(ALPHAS, function(a) data.table(
    alpha = a, observed = signif(mean(p <= a), 3), ratio = round(mean(p <= a) / a, 2)
  ))))
}

cat("=========== (a) reference = ALL events (both types) ===========\n")
for (et in c("SE", "RI")) {
  q <- nl[ET == et]
  show(loo_p(q$z, q$ID, nl$z, nl$ID), sprintf("%s values vs ALL-events reference", et))
}

cat("\n=========== (b) reference = same Event Type only (FDR family) ===========\n")
for (et in c("SE", "RI")) {
  q <- nl[ET == et]; r <- nl[ET == et]
  show(loo_p(q$z, q$ID, r$z, r$ID), sprintf("%s values vs %s-only reference", et, et))
}

cat("\n=========== (c) reference = Event Type x n-tertile ===========\n")
for (et in c("SE", "RI")) for (nt in levels(nl$n_tert)) {
  q <- nl[ET == et & n_tert == nt]; if (!nrow(q)) next
  show(loo_p(q$z, q$ID, q$z, q$ID), sprintf("%s / %s vs own-stratum reference", et, nt))
}

cat("\n=========== do the two types' z tails actually differ? ===========\n")
print(nl[, .(events = uniqueN(ID), n = .N,
             q99 = round(quantile(z, .99), 3),
             q999 = round(quantile(z, .999), 3),
             q9999 = round(quantile(z, .9999), 3)), by = ET])
cat("\nKS test SE-z vs RI-z:",
    signif(suppressWarnings(ks.test(nl[ET == "SE", z], nl[ET == "RI", z])$p.value), 3), "\n")

cat("\n=========== resolution AT FULL SCALE per FDR family ===========\n")
full <- scr[, .(events = .N), by = ET]
full[, null_z_at_200 := events * 200]
full[, floor_p := 1 / (null_z_at_200 + 1)]
full[, bh_target := c(SE = 6.64e-4, RI = 2.8e-4)[ET]]
full[, resolvable := floor_p < bh_target]
print(full)
cat("\nIf resolvable == TRUE, stratifying the reference by Event Type costs nothing:\n",
    "the family's own nulls already give more resolution than FDR requires.\n")
