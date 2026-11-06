# Makefile for phm_ensemble
#
# Covers the six set-of-six examples (01-06) plus the double-ring-of-
# cliques example -- every target below just globs by filename suffix
# (*_sim.R / *_analysis.R), so double_ring_of_cliques_sim.R and
# double_ring_of_cliques_analysis.R run alongside 01-06 automatically;
# there's no separate "07" step.
#
# make check    -- verify every R package the simulation scripts and the
#                  shared analysis library (methods_lib.R) depend on is
#                  installed; print install instructions for anything missing.
# make sim      -- run every simulation_scripts/*_sim.R script, writing its
#                  output (.rds/.txt/.png) into sim_data/ (created if needed).
# make analysis -- run every analysis_scripts/*_analysis.R script against the
#                  data in sim_data/, recording all output in $(RESULTS_FILE).
#                  The PHM ensemble's M (replicates) and rho (subsample
#                  fraction) can be overridden on the command line, e.g.:
#                    make analysis M=10 RHO=0.5
#                  Defaults (M=5, RHO=1.0) match fit_phm_ensemble()'s own
#                  defaults, and also what you get running a single
#                  analysis_scripts/*_analysis.R directly with no arguments.
#                  (double_ring_of_cliques_analysis.R also calls
#                  fit_phm_ensemble() and honors the same --M=/--rho=
#                  flags; it additionally writes
#                  double_ring_of_cliques_phm_result.pdf one directory up
#                  from analysis_scripts/, i.e. right here.)
# make clean    -- remove sim_data/ and anything in it, any stray simulated
#                  data files left directly in $(SIM_DIR), $(RESULTS_FILE), and
#                  double_ring_of_cliques_phm_result.pdf.

RSCRIPT      := Rscript
SIM_DIR      := simulation_scripts
ANALYSIS_DIR := analysis_scripts
DATA_DIR     := sim_data
RESULTS_FILE := analysis_results.txt

# PHM ensemble parameters for `make analysis`, overridable on the command
# line (e.g. `make analysis M=10 RHO=0.5`). Defaults match
# fit_phm_ensemble()'s own defaults (M=5, rho=1.0).
M   ?= 5
RHO ?= 1.0

.PHONY: check help sim analysis clean

help:
	@echo "Targets:"
	@echo "  make check     - verify required R packages are installed"
	@echo "  make sim       - generate all simulated datasets into $(DATA_DIR)/"
	@echo "  make analysis  - analyze $(DATA_DIR)/, recording output in $(RESULTS_FILE)"
	@echo "                   (override PHM ensemble params: make analysis M=10 RHO=0.5)"
	@echo "  make clean     - remove $(DATA_DIR)/ and all other generated output"

check:
	@command -v $(RSCRIPT) >/dev/null 2>&1 || { \
		echo "R does not appear to be installed  (no 'Rscript' on PATH)."; \
		echo ""; \
		echo "Install R first, then re-run 'make check':"; \
		echo "  macOS (Homebrew): brew install r"; \
		echo "  macOS (installer): https://cran.r-project.org/bin/macosx/"; \
		echo "  Linux: https://cran.r-project.org/bin/linux/"; \
		exit 1; \
	}
	@$(RSCRIPT) check_packages.R

sim:
	@mkdir -p $(DATA_DIR)
	@top="$$(pwd)"; \
	for f in $(SIM_DIR)/*_sim.R; do \
		echo "== Running $$f -> $(DATA_DIR)/ =="; \
		( cd $(DATA_DIR) && $(RSCRIPT) "$$top/$$f" ) || exit 1; \
	done; \
	echo ""; \
	echo "Done -- simulated data (.rds/.txt/.png) written to $(DATA_DIR)/"

analysis:
	@test -d $(DATA_DIR) || { \
		echo "No $(DATA_DIR)/ directory found -- run 'make sim' first."; \
		exit 1; \
	}
	@top="$$(pwd)"; \
	echo "PHM ensemble method comparison -- analysis run on $$(date) (M=$(M), rho=$(RHO))" > $(RESULTS_FILE); \
	status=0; \
	for f in $(ANALYSIS_DIR)/*_analysis.R; do \
		echo "== Running $$f (M=$(M), rho=$(RHO)) =="; \
		{ echo ""; echo "===== $$f ====="; } >> $(RESULTS_FILE); \
		if ( cd $(ANALYSIS_DIR) && $(RSCRIPT) "$$top/$$f" --M=$(M) --rho=$(RHO) ) >> $(RESULTS_FILE) 2>&1; then \
			:; \
		else \
			echo "  !! $$f exited with an error -- see $(RESULTS_FILE)"; \
			status=1; \
		fi; \
	done; \
	echo ""; \
	echo "Done -- results recorded in $(RESULTS_FILE)"; \
	exit $$status

clean:
	@rm -rf $(DATA_DIR)
	@rm -f $(SIM_DIR)/*.rds $(SIM_DIR)/*.txt $(SIM_DIR)/*.png
	@rm -f $(RESULTS_FILE)
	@rm -f double_ring_of_cliques_phm_result.pdf
	@echo "Done -- removed $(DATA_DIR)/, any simulated data files left in $(SIM_DIR)/, $(RESULTS_FILE), and double_ring_of_cliques_phm_result.pdf"
