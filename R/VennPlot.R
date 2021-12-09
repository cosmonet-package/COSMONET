
# This script contains the utility function to show venn diagram plots between low-risk (or high-risk) groups of each screening
#
## Input: groupRisk1, groupRisk2, groupRisk3 -- character giving the risk information
#         risk -- character giving the group risk 
#
## Output: Three-set venn diagram plots

VennPlot <- function(groupRisk1, groupRisk2, groupRisk3, risk){
  
  library(VennDiagram)
  # Display Venn Diagram
  venn <- display_venn(
    x = list(groupRisk1, groupRisk2, groupRisk3), 
            category.names = c(paste0(risk,"-BMD") , paste0(risk,"-DAD") , paste0(risk,"-BMD+DAD")),
                       # Numbers
                       cex = 1.5,
                       # Set names
                       cat.cex = 1.5,
                       cat.default.pos = "outer",
                       # Font
                       fontfamily = "sans",
                       cat.fontfamily = "sans"
    )
  print(venn)
}

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}