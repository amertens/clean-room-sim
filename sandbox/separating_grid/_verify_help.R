suppressMessages(library(cleanTMLE))
db <- tools::Rd_db("cleanTMLE")
rd <- db[["run_plasmode_dq_stress.Rd"]]
txt <- paste(capture.output(tools::Rd2txt(rd)), collapse=" ")
cat("installed-help mentions NA-safety:", grepl("fewer than two replicates", txt), "\n")
