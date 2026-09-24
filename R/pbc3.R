#' Mayo Clinic primary biliary cirrhosis data used as example code
#'
#' @description The dataset originates from the Mayo Clinic trial on 
#' primary biliary cirrhosis (PBC) of the liver, carried out from 1974 to 1984. 
#' It includes data from 424 PBC patients who were referred to the Mayo Clinic 
#' within this decade and met the eligibility requirements for a randomized 
#' placebo-controlled trial of D-penicillamine. However, only the initial 312 cases 
#' from the dataset were enrolled in the randomized trial. 
#' Thus, the dataset specifically pertains to these 312 patients, 
#' for whom the data is largely complete.
#'
#' @usage data(pbc2)
#' @format A data frame with 1945 observations on the following 20 variables:
#'   \describe{
#'
#'   \item{\code{id}}{patients identifier; in total there are 312 patients.}
#'
#'   \item{\code{years}}{number of years between registration and the earlier of
#'   death, transplantation, or study analysis time.}
#'
#'   \item{\code{status}}{a factor with levels \code{alive}, \code{transplanted}
#'   and \code{dead}.}
#'
#'   \item{\code{drug}}{a factor with levels \code{placebo} and
#'   \code{D-penicil}.}
#'
#'   \item{\code{age}}{at registration in years.}
#'
#'   \item{\code{sex}}{a factor with levels \code{male} and \code{female}.}
#'
#'   \item{\code{year}}{number of years between enrollment and this visit date,
#'   remaining values on the line of data refer to this visit.}
#'
#'   \item{\code{ascites}}{a factor with levels \code{No} and \code{Yes}.}
#'
#'   \item{\code{hepatomegaly}}{a factor with levels \code{No} and \code{Yes}.}
#'
#'   \item{\code{spiders}}{a factor with levels \code{No} and \code{Yes}.}
#'
#'   \item{\code{edema}}{a factor with levels \code{No edema} (i.e. no edema and
#'   no diuretic therapy for edema), \code{edema no diuretics} (i.e. edema
#'   present without diuretics, or edema resolved by diuretics), and
#'   \code{edema despite diuretics} (i.e. edema despite diuretic therapy).}
#'
#'   \item{\code{serBilir}}{serum bilirubin in mg/dl.}
#'
#'   \item{\code{serChol}}{serum cholesterol in mg/dl.}
#'
#'   \item{\code{albumin}}{albumin in mg/dl.}
#'
#'   \item{\code{alkaline}}{alkaline phosphatase in U/liter.}
#'
#'   \item{\code{SGOT}}{SGOT in U/ml.}
#'
#'   \item{\code{platelets}}{platelets per cubic ml/1000.}
#'
#'   \item{\code{prothrombin}}{prothrombin time in seconds.}
#'
#'   \item{\code{histologic}}{histologic stage of disease.}
#'
#'   \item{\code{status2}}{a numeric vector with the value 1 denoting if the
#'   patient was dead, and 0 if the patient was alive or transplanted.}
#'   
#'   \item{\code{status3}}{a numeric vector with the value 1 denoting if the
#'   patient was dead or transplanted, and 0 if the patient was alive.}
#'   
#'   \item{\code{status4}}{a numeric vector with the value 1 denoting if the
#'   patient was transplanted, and 0 if the patient was dead. Used for competing risks.}
#'   
#'   \item{\code{status5}}{a numeric vector with the value 2 if the patient was transplanted,
#'   1 denoting if the patient was dead, and 0 if the patient was alive. 
#'   Used for competing risks with censored.}
#'   
#'   \item{\code{Tyears1}}{a numeric vector with a transformed value of time-to-event outcome.}
#'
#'   \item{\code{Tyears2}}{a numeric vector with a transformed value of time-to-event outcome.}
#'   
#'   \item{\code{Tyears3}}{a numeric vector with a transformed value of time-to-event outcome.}
#'   
#'   \item{\code{Tyears4}}{a numeric vector with a transformed value of time-to-event outcome.}
#'   
#'   }
#' @keywords datasets
#' @source \code{\link[survival]{pbc}}.
#' @references
#'
#' Fleming T, Harrington D. \emph{Counting Processes and Survival Analysis}.
#' 1991; New York: Wiley.
#'
#' Therneau T, Grambsch P. \emph{Modeling Survival Data: Extending the Cox
#' Model}. 2000; New York: Springer-Verlag.
"pbc3"
