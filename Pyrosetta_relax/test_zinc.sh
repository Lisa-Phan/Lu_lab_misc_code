
#Sample constraint file for Zn binding protein
#Syntax seems like
#AtomPair Atomtype resnumberchainnumber? metal resnumberchainumber (harmonic?) (distance) (weight?)
#Angle Atomtype resnumberchainnumber? metal resnumberchainnumber Atomtype2 resnumberchainnumber (angle) (weight)

#What does this look like for heme?

AtomPair SG 5A ZN 31A HARMONIC 2.3 0.2
AtomPair SG 8A ZN 31A HARMONIC 2.3 0.2
AtomPair NE2 21A ZN 31A HARMONIC 2.0 0.2
AtomPair NE2 26A ZN 31A HARMONIC 2.0 0.2
Angle SG 5A ZN 31A SG 8A HARMONIC 1.967 0.1
Angle SG 5A ZN 31A NE2 21A HARMONIC 1.917 0.1
Angle SG 5A ZN 31A NE2 26A HARMONIC 1.844 0.1
Angle SG 8A ZN 31A NE2 21A HARMONIC 1.863 0.1
Angle SG 8A ZN 31A NE2 26A HARMONIC 1.933 0.1
Angle NE2 21A ZN 31A NE2 26A HARMONIC 1.942 0.1
