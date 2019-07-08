\name{MSanalysis}
\alias{MSanalysis}
\title{The main function of the rCPDMS package}
\description{
  This rCPDMS R-package enables automatic data analysis of multiple chemoproteomics workflows. 
}
\details{
  Current data package is able to handle the following experimental workflows:
  \enumerate{
      \item Pull-Down Assay - test the protein intensity at various dose points to identify interesting protein bindings.
      \item Cellular Thermal Shift Assay (CETSA) - test the protein intensity at various temperature points to understand the thermal stability of proteins.
      \item 2D-CETSA - test protein intensity at various dose cross various temperature to get the global picture of the protein binding information.
      % \item Pro10K - try to understand the compound-protein functional interaction.
  }
}

\section{Analysis Procedure}{
The R package will take a uniformed input data for all workflows, and perform different statistical analysis according to workflow and experimental design. Finally FDR correction was performed for p-Values given by all models. 
The input data contains the intensity related readout (intensity/LFQ/iBAQ) as well as experiment design Meta information. It is transformed from the standard Maxquant protein files or other MS protein Quantification tools (i.e. IP2) based on a mapping Meta information file. 
The data analysis includes two steps: First, data preparation and quality control (QC) and then significant protein analysis as shown below: 

  \subsection{Step 1: Data Preparation and QC}{
    \itemize{
      \item Outlier Detection: flag the outlier data points. 
      \item Normalization – Factor Based Normalization: 
          \itemize{
            \item Assume the median of each sample protein intensity is equivalent for the same or similar experiment condition.
            \item Assume protein intensity is decrease with temperature increase
          }
      \item Missing Value Imputation:
          \itemize{
            \item Assume missing value is due to below limit of detection (LOD)
            \item Use Truncated Normal Distribution to sample missing value
          }
    }
  }
 
 \subsection{Step 2: Statistical Analysis}{
    \itemize{
      \item Pull-Down
          \itemize{
            \item Linear regression \eqn{ log(Result\_Value) \sim log(Dose)}
            \item Oneway ANOVA: compare results between each dose to the specified baseline. 
            \item 4-parameter logistic curve: \eqn{ log(Result\_Value) \sim log(Dose)}, fit when more than 4 data points and at least 4 doses present. 
            \item Significant B: Draw the ratio-intensity distribution based on overall samples. 
                  Outliers are consider to be significant. 
                  Perform fitting if less than 4 doses and single repeat. 
            \item Limma: A combined approach using linear modeling and differential expression.
          }
      \item CETSA
          \itemize{
            \item 4-parameter log-logistic model: for each SampleGroup, 
                  fit \eqn{ log(Result\_Value) \sim Temperature}
          }
      \item 2D-CETSA
          \itemize{
            \item Fit 4-parameter log-logistic model to Temperature variable for each Dose (similar to CETSA analysis)
            \item Fit the same models for pull-down to Dose variable for each temperature.
          }
%      \item Pro10K
%          \itemize{
%            \item Oneway ANOVA: compare results between each SampleGroup to the specified baseline. 
%          }
    }
 }
}

\usage{
  MSanalysis(workflow, dataPath, resultType, baseline, saveFolder, 
              outlier, normalization, imputation)
}
\arguments{
  \item{workflow}{Workflow options: pulldown, CETSA, 2dCETSA}
  \item{dataPath}{File path of the input dataset}
  \item{resultType}{Intensity; LFQ; iBAQ}
  \item{baseline}{default 0 pulldown; default DMSO for other workflows}
  \item{saveFolder}{Folder path to save all outputs}
  \item{outlier}{TRUE or FALSE; whether to add outlier flag column}
  \item{normalization}{TRUE or FALSE; whether to perform normalization}
  \item{imputation}{TRUE or FALSE; whether to perform missing data imputation}
}
\value{
  A data frame of fitting results. Header may vary across different workflows. 
}
\examples{
  \dontrun{
    workflow = "CETSA" 
    dataPath = "data/CETSA_Sample.txt"
    resultType = "Intensity" 
    baseline = "DMSO"
    saveFolder = file.path("output",workflow)
    outlier = TRUE
    normalization = TRUE
    imputation = TRUE
    results = rCPDMS::MSanalysis(workflow, dataPath, resultType, baseline, saveFolder,
              outlier, normalization, imputation)
  }
}