library(data.table)
library(httr)
nonzero_cCRE_accessions <- readLines('nonzero_cCRE_accessions.txt')
tf_dt <- rbindlist(sapply(nonzero_cCRE_accessions, function(this_accession) {
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
    return(data.table(name='ERROR', n = 0L, total = 0L))
  } else {
    rbindlist(content(r, "parsed", "application/json")[[this_accession]]$tf)
  }
}), idcol = 'accession')
fwrite(tf_dt, 'tf_dt.csv')