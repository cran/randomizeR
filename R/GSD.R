#' @include OBF.R
#' @include POC.R
NULL

#' Simulate Group Sequential Design according to O'Brien and Fleming
#'
#' Calculates the mean type-I-error and go-no-go value for each stage of a group
#'  sequential design study according to O'Brien and Fleming's design.
#'
#' @param n total sample size
#' @param norep number of simulations to be conducted
#' @param K number of stages. (Currently available: \code{2,3} and \code{4})
#' @param RP the randomization procedure used.(Currently available: \code{'RAR'},
#'  \code{'PBD2'}, \code{'BSD'} - BSD(3),BSD(5) and BSD(10), \code{'EBC'} - EBC(0.5),
#'   EBC(2/3) and EBC(0.8), \code{'CHEN'} - CHEN(2/3,3), CHEN(2/3,5) and CHEN(2/3,10) )
#' @param approach the approach used, where \code{1} corresponds to allocating all patients
#'  at the begining of the trial and \code{2} corresponds to allocating patients for each
#'   stage separately (Only needed for \code{'RAR', 'BSD'} and \code{'CHEN'})
#' @examples
#' # Simulate a GSD according to O'Brien and Fleming's design with 30 patients, 10 simulation runs,
#' # 3 Stages using Random Allocation Rule as a randomization procedure and applying it to
#' # all patients at the beginning
#' GSD_OBF(30,10,2,'RAR',1)
#' @returns A Dataframe with the mean type-I-error and go-no-go values for each stage of
#' a group sequential design study according to O'Brien and Fleming's design
#' @export
#'

GSD_OBF <- function(n,norep,K,RP,approach) {
  switch (RP,
          'RAR' = switch(approach,
                         '1' = switch(K,
                                      '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                      '2' = OBF_2_RAR_A1(n,norep),
                                      '3' = OBF_3_RAR_A1(n,norep),
                                      '4' = OBF_4_RAR_A1(n,norep),
                                      "GSD with only 2,3, or 4 stages are supported at this moment"
                         ),
                         '2' = switch(K,
                                      '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                      '2' = OBF_2_RAR_A2(n,norep),
                                      '3' = OBF_3_RAR_A2(n,norep),
                                      '4' = OBF_4_RAR_A2(n,norep),
                                      "GSD with only 2,3, or 4 stages are supported at this moment"
                         )),
          'PBD2' = switch(K,
                          '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                          '2' = OBF_2_PBD2(n,norep),
                          '3' = OBF_3_PBD2(n,norep),
                          '4' = OBF_4_PBD2(n,norep),
                          "GSD with only 2,3, or 4 stages are supported at this moment"
          ),
          'BSD'=switch(approach,
                       '1' = switch(K,
                                    '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                    '2' = OBF_2_BSD_A1(n,norep),
                                    '3' = OBF_3_BSD_A1(n,norep),
                                    '4' = OBF_4_BSD_A1(n,norep),
                                    "GSD with only 2,3, or 4 stages are supported at this moment"
                       ),
                       '2' = switch(K,
                                    '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                    '2' = OBF_2_BSD_A2(n,norep),
                                    '3' = OBF_3_BSD_A2(n,norep),
                                    '4' = OBF_4_BSD_A2(n,norep),
                                    "GSD with only 2,3, or 4 stages are supported at this moment"
                       )

          ),
          'EBC' = switch(K,
                         '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                         '2' = OBF_2_EBC(n,norep),
                         '3' = OBF_3_EBC(n,norep),
                         '4' = OBF_4_EBC(n,norep),
                         "GSD with only 2,3, or 4 stages are supported at this moment"
          ),
          'CHEN'=switch(approach,
                        '1' = switch(K,
                                     '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                     '2' = OBF_2_CHEN_A1(n,norep),
                                     '3' = OBF_3_CHEN_A1(n,norep),
                                     '4' = OBF_4_CHEN_A1(n,norep),
                                     "GSD with only 2,3, or 4 stages are supported at this moment"
                        ),
                        '2' = switch(K,
                                     '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                     '2' = OBF_2_CHEN_A2(n,norep),
                                     '3' = OBF_3_CHEN_A2(n,norep),
                                     '4' = OBF_4_CHEN_A2(n,norep),
                                     "GSD with only 2,3, or 4 stages are supported at this moment"
                        )

          ),
          "GSD with just RAR, CHEN, EBC, BSD and PBD(2) is supported at this moment"
  )
}

#' Simulate Group Sequential Design according to Pocock
#'
#' Calculates the mean type-I-error and go-no-go value for each stage of a group
#'  sequential design study according to Pocock's design.
#'
#' @param n total sample size
#' @param norep number of simulations to be conducted
#' @param K number of stages. (Currently available: \code{2,3} and \code{4})
#' @param RP the randomization procedure used.(Currently available: \code{'RAR'},
#'  \code{'PBD2'}, \code{'BSD'} - BSD(3),BSD(5) and BSD(10), \code{'EBC'} - EBC(0.5),
#'   EBC(2/3) and EBC(0.8), \code{'CHEN'} - CHEN(2/3,3), CHEN(2/3,5) and CHEN(2/3,10) )
#' @param approach the approach used, where \code{1} corresponds to allocating all patients
#'  at the begining of the trial and \code{2} corresponds to allocating patients for each
#'   stage separately (Only needed for \code{'RAR', 'BSD'} and \code{'CHEN'})
#' @examples
#' # Simulate a GSD according to Pocock's design with 30 patients, 10 simulation runs,
#' # 2 Stages using Random Allocation Rule as a randomization procedure and applying it to
#' # all patients at the beginning
#' GSD_POC(30,10,2,'RAR',1)
#'
#' @returns A Dataframe with the mean type-I-error and go-no-go values for each stage of
#' a group sequential design study according to Pocock's design/
#' @export
#'

GSD_POC <- function(n,norep,K,RP,approach) {
  switch (RP,
          'RAR' = switch(approach,
                         '1' = switch(K,
                                      '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                      '2' = POC_2_RAR_A1(n,norep),
                                      '3' = POC_3_RAR_A1(n,norep),
                                      '4' = POC_4_RAR_A1(n,norep),
                                      "GSD with only 2,3, or 4 stages are supported at this moment"
                         ),
                         '2' = switch(K,
                                      '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                      '2' = POC_2_RAR_A2(n,norep),
                                      '3' = POC_3_RAR_A2(n,norep),
                                      '4' = POC_4_RAR_A2(n,norep),
                                      "GSD with only 2,3, or 4 stages are supported at this moment"
                         )
          ),
          'PBD2' = switch(K,
                          '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                          '2' = POC_2_PBD2(n,norep),
                          '3' = POC_3_PBD2(n,norep),
                          '4' = POC_4_PBD2(n,norep),
                          "GSD with only 2,3, or 4 stages are supported at this moment"
          ),
          'BSD'=switch(approach,
                       '1' = switch(K,
                                    '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                    '2' = POC_2_BSD_A1(n,norep),
                                    '3' = POC_3_BSD_A1(n,norep),
                                    '4' = POC_4_BSD_A1(n,norep),
                                    "GSD with only 2,3, or 4 stages are supported at this moment"
                       ),
                       '2' = switch(K,
                                    '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                    '2' = POC_2_BSD_A2(n,norep),
                                    '3' = POC_3_BSD_A2(n,norep),
                                    '4' = POC_4_BSD_A2(n,norep),
                                    "GSD with only 2,3, or 4 stages are supported at this moment"
                       )

          ),
          'EBC' = switch(K,
                         '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                         '2' = POC_2_EBC(n,norep),
                         '3' = POC_3_EBC(n,norep),
                         '4' = POC_4_EBC(n,norep),
                         "GSD with only 2,3, or 4 stages are supported at this moment"
          ),
          'CHEN'=switch(approach,
                        '1' = switch(K,
                                     '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                     '2' = POC_2_CHEN_A1(n,norep),
                                     '3' = POC_3_CHEN_A1(n,norep),
                                     '4' = POC_4_CHEN_A1(n,norep),
                                     "GSD with only 2,3, or 4 stages are supported at this moment"
                        ),
                        '2' = switch(K,
                                     '1' = print("GSD with only 2,3, or 4 stages are supported at this moment"),
                                     '2' = POC_2_CHEN_A2(n,norep),
                                     '3' = POC_3_CHEN_A2(n,norep),
                                     '4' = POC_4_CHEN_A2(n,norep),
                                     "GSD with only 2,3, or 4 stages are supported at this moment"
                        )

          ),
          "GSD with just RAR, CHEN, EBC, BSD and PBD(2) is supported at this moment"
  )
}
