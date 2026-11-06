for K in 1 2 3 4 5 6;
do
  echo "Running admixture for K$K..."
  admixture --cv /labgenomaarea2/saguis-lucero/2_DataExploration/LW_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac.ped $K | tee log_LW_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac.${K}.out
done

        #check the CV values over all the log files from the command above
grep -h CV log_LW_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac*.out
for K in 1 2 3 4 5 6;
do
  echo "Running admixture for K$K..."
  admixture --cv /labgenomaarea2/saguis-lucero/2_DataExploration/TI_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac.ped $K | tee log_TI_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac.${K}.out
done

        #check the CV values over all the log files from the command above
grep -h CV log_TI_1SNPevery10kb_hwe0.001maf0.05_Missing30OnlyBialelicsQuali30_minDP10maxDP40_OnlyAut_NormRawQq30d100_refcalJac*.out
