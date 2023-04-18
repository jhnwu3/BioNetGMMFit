python3 graph.py -f test/l6p/l6p_estimates.txt -g CI_truth -n "Linear 6 Protein" -r test/l6p/tr.csv 
python3 graph.py -f test/l6p/l6p_est_m.csv -g CI_truth -n "Linear 6 Protein Means Only" -r test/l6p/tr.csv 
python3 graph.py -f test/l6p/l6p_est_mv.csv -g CI_truth -n "Linear 6 Protein Means and Variances Only" -r test/l6p/tr.csv 
python3 graph.py -f test/l6p/6_pro_lint1.50_leastCostMoments.csv -g Moments -n "Linear 6 Protein" -m 6
python3 graph.py -f test/l6p/l6p_mv.csv -g Moments -n "Means and Variances Only" -m 6
python3 graph.py -f test/l6p/l6p_m.csv -g Moments -n "Means Only" -m 6

python3 graph.py -f test/yeast/yeast_estimates.csv -g CI_truth -n "Yeast" -r test/yeast/true_rates.csv
python3 graph.py -f test/yeast/yeastt1.00_leastCostMoments.csv -g Moments -n "Yeast at t=1" -m 5
python3 graph.py -f test/yeast/yeastt5.00_leastCostMoments.csv -g Moments -n "Yeast at t=5" -m 5