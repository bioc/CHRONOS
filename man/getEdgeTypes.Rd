\name{getEdgeTypes}
\alias{getEdgeTypes}
\title{Map various types of gene-gene interactions in KGML files to edge types 
in corresponding pathway graphs.}
\usage{
getEdgeTypes(type)
}
\arguments{
\item{type}{A vector of interaction types.}
}
\value{
If an interaction type has been supplied, the corresponging edge 
types are returned. 
If not, the complete mapping is returned.
}
\description{
Map various types of gene-gene interactions in KGML files to edge types 
in corresponding pathway graphs.
}
\details{
Edge types

\tabular{llllllll}{
    activation \tab 1 \tab inhibition \tab 2 \tab apathetic \tab 3 \tab 
    no interaction \tab 4
} 

Default interaction - edge type mapping

\tabular{rrrrrr}{
    01 \tab unknown \tab 3 \tab 02 \tab activation \tab 1\cr
    03 \tab inhibition \tab 2 \tab 04 \tab binding/association \tab 3\cr
    05 \tab expression \tab 1 \tab 06 \tab repression \tab 2\cr
    07 \tab phosphorylation \tab 3 \tab 08 \tab dephosphorylation \tab 3\cr
    09 \tab ubiquitination \tab 3 \tab 10 \tab dissociation \tab 3\cr
    11 \tab indirect effect \tab 3 \tab 12 \tab state change \tab 3\cr
    13 \tab compound \tab 3 \tab 14 \tab hidden compound \tab 3\cr
    16 \tab missing interaction 
        \tab 3 \tab 16 \tab activation_phosphorylation \tab 1\cr
    17 \tab activation_dephosphorylation 
        \tab 1 \tab 18 \tab activation_ubiquitination \tab 1\cr
    19 \tab activation_indirect effect 
        \tab 1 \tab 20 \tab activation_binding/association \tab 1\cr
    21 \tab activation_inhibition 
        \tab 3 \tab 22 \tab activation_methylation \tab 1\cr
    23 \tab inhibition_phosphorylation 
        \tab 2 \tab 24 \tab inhibition_dephosphorylation \tab 2\cr
    25 \tab inhibition_ubiquitination 
        \tab 2 \tab 26 \tab inhibition_indirect effect \tab 2\cr
    27 \tab inhibition_binding/association 
        \tab 2 \tab 28 \tab inhibition_expression \tab 2\cr
    29 \tab inhibition_methylation 
        \tab 2 \tab 30 \tab compound_expression \tab 1\cr
    31 \tab compound_activation 
        \tab 1 \tab 32 \tab compound_inhibition \tab 2\cr
}

\tabular{rrrrrr}{
    33 \tab compound_activation_indirect effect \tab 1 \cr 
    34 \tab compound_activation_phosphorylation \tab 1\cr
    35 \tab phosphorylation_indirect effect \tab 3 \cr 
    36 \tab phosphorylation_binding/association \tab 3\cr
    37 \tab phosphorylation_dissociation \tab 3 \cr 
    38 \tab dephosphorylation_indirect effect \tab 3\cr
    39 \tab binding/association_missing interaction \tab 3 \cr 
    40 \tab binding/association_indirect effect \tab 3\cr
    41 \tab expression_indirect effect \tab 1 \cr 
    42 \tab repression_indirect effect \tab 2\cr
    43 \tab ubiquitination_inhibition \tab 2 \cr 
    44 \tab dissociation_missing interaction \tab 3\cr
    45 \tab indirect effect_phosphorylation \tab 3  \tab \tab \cr
    46 \tab activation_phosphorylation_binding/association \tab 1 \cr
    47 \tab activation_phosphorylation_indirect effect \tab 1\cr
}
}
\examples{

# Example 1

# Retreive edge types for phosphorylation and dephosphorylation.
getEdgeTypes(c(7,8)) 

# Example 2

# Returns a data frame containing the innteraction - edge type mapper.
types <- getEdgeTypes()

# Set phosphorylation to inhibition.
types[8,2] <- 2


}

