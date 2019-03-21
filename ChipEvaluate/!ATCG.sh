awk 'BEGIN{FS="\t"} {if ($NF=="1"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="2"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="3"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="4"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="5"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="6"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="7"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="8"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt
awk 'BEGIN{FS="\t"} {if ($NF=="9"){print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 | grep -v A/T|grep -v C/G |wc -l >>notATCG.txt

