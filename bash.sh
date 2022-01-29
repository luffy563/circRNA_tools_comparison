# # find_circ
{ time nohup bash bash_FC.sh &>output/logs/FC.log & } 2>output/logs/FC.time

# # CIRCexplorer2
{ time nohup bash bash_CIRC2-1.sh &>output/logs/CIRC2.log & } 2>output/logs/CIRC2.time

# # CircRNAfinder
{ time nohup bash bash_CF.sh &>output/logs/CF.log & } 2>output/logs/CF.time

# # CDBG
{ time nohup bash bash_CDBG.sh &>output/logs/CDBG.log & } 2>output/logs/CDBG.time

# # CircMarker
{ time nohup bash bash_CM.sh &>output/logs/CM.log & } 2>output/logs/CM.time

# # CIRI2
{ time nohup bash bash_CIRI2.sh &>output/logs/CIRI2.log & } 2>output/logs/CIRI2.time

# # DCC
{ time nohup bash bash_DCC.sh &>output/logs/DCC.log & } 2>output/logs/DCC.time

# # Mapsplice
{ time nohup bash bash_MP.sh &>output/logs/MP.log & } 2>output/logs/MP.time

# # Segemeghl
{ time nohup bash bash_SE.sh &>output/logs/SE.log & } 2>output/logs/SE.time

# # KNIFE
{ time nohup bash bash_KF.sh &>output/logs/NCL.log & } 2>output/logs/NCL.time

# # UROBORUS
{ time nohup bash bash_UB.sh &>output/logs/UB.log & } 2>output/logs/UB.time

