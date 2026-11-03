for K in 1 2 3 4 5 6;
do
  echo "Running admixture for K$K..."
  admixture --cv /labgenomaarea2/saguis-lucero/2_DataExploration/LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.ped $K | tee log_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.${K}.out
done

        #check the CV values over all the log files from the command above
grep -h CV log_LW_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac*.out
for K in 1 2 3 4 5 6;
do
  echo "Running admixture for K$K..."
  admixture --cv /labgenomaarea2/saguis-lucero/2_DataExploration/TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.ped $K | tee log_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac.${K}.out
done

        #check the CV values over all the log files from the command above
grep -h CV log_TI_HWE001_1SNPevery10kb_minDP13_OnlyAut_NoMissMAF0.05Norm_MinDepthFilter_BiallelicQq30_Saguis_refcalJac*.out
