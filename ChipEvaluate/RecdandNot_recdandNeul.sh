awk 'BEGIN{FS=="\t"}{if ($NF=="9" && $(NF-2)!="not_recommended" && $(NF-2)!="neutral") {print $0}}' Axiom_KPSmilet_redo_scored.txt.sorted3600 > 3600_9.txt
