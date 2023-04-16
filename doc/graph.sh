python3 graph.py -f test/l6p/l6p_estimates.txt -g CI_truth -n "Linear 6 Protein" -r test/l6p/tr.csv
python3 graph.py -f test/l6p/6_pro_lint1.50_leastCostMoments.csv -g Moments -n "Linear 6 Protein" 

python3 graph.py -f test/yeast/yeast_estimates.csv -g CI_truth -n "Yeast" -r test/yeast/true_rates.csv
python3 graph.py -f test/l6p/6_pro_lint1.50_leastCostMoments.csv -g Moments -n "Linear 6 Protein" 