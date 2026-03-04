.PHONY: run run-quick knit test clean

# Run full pipeline (B=500 bootstrap)
run:
	Rscript scripts/00_run_all_stages.R

# Run pipeline with reduced bootstrap (B=200)
run-quick:
	Rscript scripts/00_run_all_stages.R --quick

# Knit the report (requires outputs/ from 'make run')
knit:
	cd analysis && quarto render case_study_clean_room.qmd --to html

# Run smoke test
test:
	Rscript scripts/99_smoke_test.R

# Remove generated outputs (preserves source code)
clean:
	rm -rf outputs/stage1 outputs/stage2 outputs/stage3 outputs/report
	rm -f outputs/*.json outputs/*.csv outputs/*.rds
	rm -f logs/*.txt logs/*.log
