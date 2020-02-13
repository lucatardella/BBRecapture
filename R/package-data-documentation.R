#' Great-Copper butterfly data
#'
#' Great Copper butterflies in Willamette Valley of Oregon. Once presumed to be extinct from western Oregon since 1970 (R. A. Pyle, 2002, Butterflies of Cascadia, Seattle Audubon Society), the great copper (Lycaena xanthoides Lycaenidae), was rediscovered in several wetland prairie remnants in the summer of 2004 (see Severns and Villegas 2005 for an account). Until this rediscovery, the habitat and host plants were unknown for western Oregon populations, and the great copper was conspicuously sparse in collections, with a total of twelve specimens captured between 1896 and 1970.  A current estimate of the extant populations is alarmingly low, about 100 individuals across three subpopulations in the Willamette Valley. The great copper occurs in much larger numbers through much of California. \url{http://people.oregonstate.edu/~wilsomar/Persp_GrtCop.htm}
#'
#'This dataset has been first analyzed in Ramsey F, Severns P (2010) and later re-analyzed in Farcomeni (2011) and Alunni Fegatelli and Tardella (2012)
#'
#' @docType data
#'
#' @usage data(greatcopper)
#'
#' @format A matrix with 45 rows (observed butterflies) and 8 columns (capture occasions)
#'
#' @keywords Datasets
#'
#' @references
#'
#' Severns, P.M. and Villegas, S. (2005) Butterflies Hanging on to Existence in the Willamette Valley: A Relict Population of the Great Copper (Lycaena xanthoides Boisduval), Northwest Science, Vol. 79, No. 1, 77--80.
#'
#' Severns, P.M. and Villegas, S. (2006) Conserving a wetland butterfly: quantifying early lifestage survival through seasonal flooding, adult nectar, and habitat preference, Journal of Insect Conservation 10 (4), 361-370
#'
#' Ramsey F, Severns P (2010) Persistence models for mark-recapture. Environmental and Ecological Statistics 17(1):97--109
#'
#' Farcomeni A. (2011) Recapture models under equality constraints for the conditional capture probabilities. Biometrika 98(1):237--242
#'
#' Alunni Fegatelli, D. and Tardella, L. (2012) Improved inference on capture recapture models with behavioural effects. Statistical Methods & Applications Volume 22, Issue 1, pp 45-66 10.1007/s10260-012-0221-4
#'
#'
#' @source Ramsey F, Severns P (2010) Persistence models for mark-recapture. Environmental and Ecological Statistics 17(1):97--109
#'
#' @examples
#' data(greatcopper)
#' str(greatcopper)
#'
"greatcopper"
NULL
#' Mouse (Microtus Pennsylvanicus) Dataset
#'
#' The mouse (Microtus pennsylvanicus) data were first discussed in Nichols, Pollock and Hines (1984). The original live-trapping experiment was conducted monthly from June to December, 1980. During each month, the capture-recapture procedure was repeated for 5 consecutive days.
#' The detailed data are given in Williams, Nichols and Conroy (2002, pp. 525-528). We use the data collected in June. A total of 104 distinct mice were caught in the experiment.
#'
#' @docType data
#'
#' @usage data(mouse)
#'
#' @format A matrix with 104 rows (observed animals) and 5 columns (capture occasions)
#'
#' @keywords Datasets
#'
#' @references
#'
#' J. D. Nichols, K. H. Pollock and J. E. Hines, The Use of a Robust Capture-Recapture Design in Small Mammal Population Studies: A Field Example with Microtus Pennsylvanicus, Acta Theriologica, vol. 29. 30:357-365, 1984
#'
#' @source J. D. Nichols, K. H. Pollock and J. E. Hines, The Use of a Robust Capture-Recapture Design in Small Mammal Population Studies: A Field Example with Microtus Pennsylvanicus, Acta Theriologica, vol. 29. 30:357-365, 1984
#'
#' @examples
#' data(mouse)
#' str(mouse)
#'
"mouse"
NULL
#' Flat-tailed Horned Lizard Dataset
#'
#' Data from multiple searches for flat-tailed horned lizards (\emph{Phrynosoma mcalli}) on a plot in Arizona, USA.
#'
#' @details
#'
#' The flat-tailed horned lizard (\emph{Phrynosoma mcalli}) is a desert lizard found in parts of southwestern Arizona, southeastern California and northern Mexico. There is considerable concern about its conservation status. The species is cryptically colored and has the habit of burying under the sand when approached, making it difficult or impossible to obtain a complete count (Grant and Doherty 2007).
#'
#' A total of 68 individuals were captured 134 times. Exactly half of the individuals were recaptured exactly only once. This dataset is also includedin the \code{secr} package.
#'
#' @docType data
#'
#' @usage data(hornedlizard)
#'
#' @format A matrix with 68 rows (observed lizards) and 14 columns (capture occasions)
#'
#' @keywords Datasets
#'
#' @references
#'
#' Grant, T. J. and Doherty, P. F. (2007) Monitoring of the flat-tailed horned lizard with methods incorporating detection probability. \emph{Journal of Wildlife Management} \bold{71}, 1050--1056
#'
#' Royle, J. A. and Young, K. V. (2008) A hierarchical model for spatial capture–recapture data. Ecology 89, 2281–2289.
#'
#' @source Royle, J. A. and Young, K. V. (2008)
#'
#' @examples
#' data(hornedlizard)
#' str(hornedlizard)
#'
"hornedlizard"
NULL
#' Lizard Dataset
#'
#' Giant day gecko's capture history data observed in Masoala rainforest exhibit, Zurich Zoo. This datatset was first ussed in Wanger T.C. et al. (2009).
#'
#' @details
#'
#' The giant day gecko (Phelsuma madagascariensis grandis) is a tropical reptile naturally occurring in the northern part of Madagascar. Data represent an extensive capture-recapture experiment at the  Masoala rainforest exhibit, Zurich Zoo. The researchers used the individual color patterns of this gecko species for photo recognition. Due to the high number of sampling sessions, this study is based on an unusually good dataset.  Moreover, it is a captive population and, therefore, we can be sure that the closed population assumption is valid. However, the dataset is comparable to natural conditions because of the dimensions of the Masoala rainforest exhibit; it is currently the second largest tropical exhibit in the world. Additional variation in the dataset results from time-dependent variation in recapture rates (obvious dependence of gecko activity on daily weather) and the photographic capture method (it is harder to spot and photograph juvenile compared to adult geckos).
#'
#' @docType data
#'
#' @usage data(lizard)
#'
#' @format A matrix with 68 rows (observed units) and 30 columns (capture occasions)
#'
#' @keywords Datasets
#'
#' @references
#'
#' Wanger T.C. et al. (2009) How to monitor elusive lizards: comparison of capture-recapture methods on giant day geckos (Gekkonidae, Phelsuma madagascariensis grandis) in the Masoala rainforest exhibit, Zurich Zoo. Ecological Research 24(2):345--353.
#'
#' @source Wanger T.C. et al. (2009)
#'
#' @examples
#'
#' data(lizard)
#' str(lizard)
#'
"lizard"
NULL
