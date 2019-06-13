#!/bin/sed -f

# Delete all lines that are not ATOM HETATM or TER
/^ATOM\|^HETATM\|^TER/!d

# Translate ions to names Leap likes
s/^HETATM \(....\) CA    CA/ATOM   \1  CAL CAL/
s/^HETATM \(....\) MG    MG/ATOM   \1  MAG MAG/
s/^HETATM \(....\) ZN    ZN/ATOM   \1  ZIN ZIN/
s/^HETATM \(....\) CL    CL/ATOM   \1  CHL CHL/
s/^HETATM \(....\)  K     K/ATOM   \1  POT POT/

s/^HETATM \(..........\)HEM/ATOM   \1HEM/ 

# Delete all waters
# /HETATM ....  O   HOH /d

# Delete all remaining HETATM records (usually ligand+waters)
/HETATM/d
