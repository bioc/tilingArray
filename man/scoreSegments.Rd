## This algorithm scores segments.
## For each segment, we calculate and/or record its:
## * chromosome, strand, start position, end position (in bases)
## * level (mean-background)
## * p-value (one-sample t-test againt hypothesis mean=0)
##
## If it overlaps with an annotated feature
## * the name of that feature
## * the signed distance between its start and that of the feature
## * the signed distance between its end and that of the feature
## If not:
## * the distance between its start and the end of the next feature to the left
## * the distance between its end and the start of the next feature to the left
## 
## and the same information again for the opposite strand
##--------------------------------------------------------------------------------

## For the definition of pseudogenes at SGD, see Docs/PseudogenesAtSGD.pdf
