library(data.table)
library(httr)
library(pbapply)

pboptions(type='timer')

cCRE_accessions <- readLines('all_cCRE_accessions.txt')
curr_tf_dt <- fread('tf_dt.csv')
cCRE_accessions <- curr_tf_dt[name != 'ERROR', cCRE_accessions[!cCRE_accessions %in% accession]]
tf_dt <- rbindlist(pbapply::pbsapply(cCRE_accessions, function(this_accession) {
  r <-
    RETRY(verb = 'POST',
          url = 'https://screen-beta-api.wenglab.org/dataws/re_detail/tfIntersection',
          body = sprintf(
            '{"assembly": "GRCh38", "accession": "%s"}',
            this_accession
          ),
          content_type_json(),
          pause_min = 301
    )
  if (http_error(r)){
    
    print(this_accession)
    warn_for_status(r)
    res <- data.table(name='ERROR', n = 0L, total = 0L)
  } else {
    res <- rbindlist(content(r, "parsed", "application/json")[[this_accession]]$tf)
  }
  if(nrow(res) == 0) {
    res <- data.table(name='EMPTY', n = 0L, total = 0L)
  }
  res[, accession:=this_accession]
  fwrite(res[, .(accession, name, n, total)], 'tf_dt.csv', append = TRUE)
}), idcol = 'accession')
