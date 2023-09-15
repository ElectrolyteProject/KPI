# For_DFT_MD
整理日期：2022/10/05
整理人：高宇辰

本数据集记录所有进行DFT和MD的初始电解液分子
目前由于生成规则共有三套，即
1. 初始原子为C、O、F，数据集有重原子C、O、F
2. 初始原子为O，数据集有重原子C、O、F
3. 初始原子为O，不含C=O

分子数据库生成顺序：
1. Strat_from_COF：list_with_nodes1.dat - list_with_nodes7.dat
2. 2_Electrolyte_Project_StartfromO：Nodes8_exceptOH_includeO.dat
3. 4_20220619_Ether Molecular Database (non-fluorinated)：list_with_nodes9_ether_non-fluorinated_exceptOH.dat
4. 4_Nodes9_exceptOH_includeO_compared_with_ether_non-fluorinated_exceptOH：Nodes9_exceptOH_includeO_compared_with_ether_non-fluorinated_exceptOH.smi
5. 5_list_with_nodes10_ether_non-fluorinated_exceptOH_edges9：list_with_nodes10_edges9_ether_non-fluorinated_exceptOH.dat