\name{rankCells}
\alias{rankCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rank cell types under a scRNA-Seq framework
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankCells(disease,path,annotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disease}{
  %%     ~~Describe \code{x} here~~
  \code{The name of the disease under study}
  }
  \item{path}{
  \code{The path to the directory containing the data}
  }
  \item{annotate}{
  \code{Define whether the data is annotated or if annoation is required}
  }
  \item{scenario}{
  \code{Define which scenario/approcah to use (MalaCards=default, or Hypothesis-driven)}
  }
  \item{checkdrug}{
  \code{Define whether to use drug names (MalaCards approach=default), or dug mode of actions (MOAs) (Hypothesis-driven)}
  }
  \item{userlabel}{
  \code{Define the name of the SEURAT object column containing disease-control labels (only for annotated data)}
  }
  \item{usercelltype}{
  \code{Define the name of the SEURAT object column containing cell cluster IDs (only for annotated data)}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{disease }{A Character string with the name of disease}
\item{path }{A Character string with a valid path to the directory containing the scRNA-Seq data.}
\item{annotate }{A Boolean (TRUE/FALSE) denoting whether to annotate the data.}
\item{scenario }{A  Character string with the name of scenarion (default="Malacards") alternative option="Hypothesis".}
\item{checkdrug }{A Boolean (TRUE/FALSE) denoting whether use drug names (defauls=TRUE) alternative option for MOAs=FALSE.}
\item{userlabel }{A  Character string with the name of the SEURAT object column which contains the information labels with respect to the conditions of the experiment (disease-control labels).}
\item{usercelltype }{A  Character string with the name of the SEURAT object column which contains the annotated cell cluster IDs.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
