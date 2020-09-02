import os

os.system("rsync -av -e ssh --include=\"expression_pattern_best.tsv\" --include=\"gene_best.yml\" --include=\"rmse_data.tsv\" \
	--include=\"/paper_data*\" --include=\"/final*\" --exclude=\"expression_pattern_*\" --exclude=\"gene_*\" --exclude=\"config.yml\" \
	sbs2756@wilkcomp02.ccbb.utexas.edu:/stor/work/Wilke/sshah/pinetree-toys/results/ /home/sahil/pinetree-toys/results/")